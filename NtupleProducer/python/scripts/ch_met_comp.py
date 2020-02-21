import os, re
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(False)
ROOT.gStyle.SetErrorX(0.5)
ROOT.gErrorIgnoreLevel = ROOT.kWarning

from array import array
from math import pow, sin, cos, hypot, atan2, sqrt
import itertools
from FastPUPPI.NtupleProducer.scripts.makeJecs import _progress

def calcHT(jets):
    return sum(j[0] for j in jets) if jets else 0
def calcMHTPhi(jets):
    if jets:
        X = sum(j[0]*sin(j[2]) for j in jets)
        Y = sum(j[0]*cos(j[2]) for j in jets)
        return hypot(X,Y), atan2(-Y,-X) # negative vector sum
    return 0,0
# soft term is MET+jets (MET=-VectorSum(jets) + ST)
def calcST(met, jets):
    X = met[0]*sin(met[1])
    Y = met[0]*cos(met[1])
    if jets:
        X += sum(j[0]*sin(j[2]) for j in jets)
        Y += sum(j[0]*cos(j[2]) for j in jets)
    return hypot(X,Y), atan2(Y,X)
def calcMHT(jets):
    return hypot(sum(j[0]*sin(j[2]) for j in jets),sum(j[0]*cos(j[2]) for j in jets))  if jets else 0
def vectorSum(jets):
    return hypot(sum(j[0]*sin(j[1]) for j in jets),sum(j[0]*cos(j[1]) for j in jets))  if jets else 0
def calcMHTSTv1(jets, jetsUnCorr, met):
    mht = calcMHTPhi(jets)
    softTerm = calcST(met,jetsUnCorr)
    return vectorSum( [mht, softTerm] )
def calcMHTSTv2(met, corrs):
    mht_corrs = calcMHTPhi(corrs) # vector sum of corrs
    return vectorSum( [mht_corrs, met] )

def jetRes(pt,eta):
    #fast approx (in GeV), likely bad approx below ~15 GeV
    return pt * (0.08 + 1.5/pt)
def calcHTSig(jets, quantile):
    mht = calcMHT(jets)
    sigma = sqrt(sum([jetRes(j[0],j[1]) for j in jets]))
    return mht + quantile * sigma
def calcHTSigST1(jets,jetsUnCorr,met, quantile):
    mhtST = calcMHTSTv1(jets, jetsUnCorr, met)
    sigma = sqrt(sum([jetRes(j[0],j[1]) for j in jets]))
    return mhtST + quantile * sigma
def calcHTSigST2(jets, met,corrs, quantile):
    mhtST = calcMHTSTv2(met,corrs)
    sigma = sqrt(sum([jetRes(j[0],j[1]) for j in jets]))
    return mhtST + quantile * sigma

def calcMJJ(jets):
    if len(jets) <= 1: return 0
    lvclass = ROOT.Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<double>")
    jp4s = [ lvclass(j[0],j[1],j[2],0) for j in jets ]
    return max((j1+j2).M() for (j1,j2) in itertools.combinations(jp4s,2))

class CalcJ:
    def __init__(self,index):
        self._index = index
    def __call__(self,jets):
        if len(jets) <= self._index: return 0
        jsort = sorted(jets, key = lambda j : -j[0])
        return jsort[self._index][0]
class CalcJ2_MJJcut:
    def __init__(self,mjj):
        self._mjj = mjj
    def __call__(self,jets):
        if len(jets) <= 1: return 0
        lvclass = ROOT.Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<double>")
        jp4s = [ lvclass(j[0],j[1],j[2],0) for j in jets ]
        jjs = [ jj for jj in itertools.combinations(jp4s,2) if (jj[0]+jj[1]).M() > self._mjj ]
        if len(jjs) <= 0: return 0
        return max(min(jj[0].Pt(), jj[1].Pt()) for jj in jjs)

def makeCalc(what):
    if what == "ht": 
        return calcHT
    if what == "mht": 
        return calcMHT
    if what == "mjj": 
        return calcMJJ
    if re.match(r"jet\d+$", options.var): 
        return CalcJ(int(options.var.replace("jet",""))-1)
    if what.startswith("ptj-mjj"): 
        return CalcJ2_MJJcut(float(what.replace("ptj-mjj","")))

def makeGenArray(tree, what, ptCut, etaCut, _cache={}):
    _key = (id(tree),what,int(ptCut*100),int(etaCut*1000))
    if _key in _cache: return _cache[_key]
    if what == "metmht":
        met = makeGenArray(tree, "met",     0,    5.0, _cache=_cache)
        mht = makeGenArray(tree, "mht", ptCut, etaCut, _cache=_cache)
        ret = map(min, zip(met,mht))
        _cache[_key] = ret
        return ret
    if "met" in what:
        ret = makeGenMETArray(tree, what, etaCut)
        _cache[_key] = ret
        return ret
    progress = _progress("Reading GenJets ...")
    calc = makeCalc(what)
    ret = []
    tree.SetBranchStatus("*",0);
    tree.SetBranchStatus("nGenJets",1);
    tree.SetBranchStatus("GenJets_*",1);
    for i in xrange(tree.GetEntries()):
        tree.GetEntry(i)
        pt,eta,phi = tree.GenJets_pt, tree.GenJets_eta, tree.GenJets_phi
        jets = [ (pt[j],eta[j],phi[j]) for j in xrange(tree.nGenJets) if pt[j] > ptCut and abs(eta[j]) < etaCut ]
        ret.append(calc(jets))
    _cache[_key] = ret
    progress.done("done, %d entries" % len(ret))
    return ret
def makeGenMETArray(tree, what, etaCut):
    if   etaCut <= 1.5: post = "MetBarrel_pt" 
    elif etaCut <= 2.4: post = "MetCentral_pt"
    else:               post = "Met_pt"
    progress = _progress("Reading gen"+post+" ...")
    tree.SetBranchStatus("*",0);
    tree.SetBranchStatus("gen"+post,1);
    ret = []
    for i in xrange(tree.GetEntries()):
        tree.GetEntry(i)
        ret.append(getattr(tree,"gen"+post))
    progress.done("done, %d entries" % len(ret))
    return ret
def makeCorrArray(tree, what, whatOpts, obj, ptCorrCut, etaCut, corr, _cache={}):
    _key = (id(tree),what,obj,int(ptCorrCut*100),int(etaCut*1000))
    if _key in _cache: return _cache[_key]
    if "met" in what:
        ret = makeRecoMETArray(tree, what, obj, etaCut)
        _cache[_key] = ret
        return ret
    
    ret = []
    if not tree.GetBranch("n"+obj+"Jets"): 
        return None
    progress = _progress("Reading "+obj+"Jets ...")
    tree.SetBranchStatus("*",0);
    tree.SetBranchStatus("n"+obj+"Jets",1);
    tree.SetBranchStatus(obj+"Jets_pt",1);
    tree.SetBranchStatus(obj+"Jets_eta",1);
    tree.SetBranchStatus(obj+"Jets_phi",1);
    if   etaCut <= 1.5: post = "MetBarrel_pt" 
    elif etaCut <= 2.4: post = "MetCentral_pt"
    else:               post = "Met_pt"
    tree.SetBranchStatus(obj+post,1);
    tree.SetBranchStatus(obj+post[:-2]+"phi",1);

    bias=0
    if what=="met": calc = lambda jets,jetsUnCorr,corrs,met : met[0]
    elif what=="mht": calc = lambda jets,jetsUnCorr,corrs,met : calcHT(jets)
    elif what=="mhtNoJEC": calc = lambda jets,jetsUnCorr,corrs,met : calcHT(jetsUnCorr)
    elif what=="mhtST1": calc = lambda jets,jetsUnCorr,corrs,met : calcMHTSTv1(jets,jetsUnCorr,met)
    elif what=="mhtST2": calc = lambda jets,jetsUnCorr,corrs,met: calcMHTSTv2(met,corrs)
    elif what=="mhtSig_p1": calc = lambda jets,jetsUnCorr,corrs,met: calcHTSig(jets, 1)
    elif what=="mhtSig_p2": calc = lambda jets,jetsUnCorr,corrs,met: calcHTSig(jets, 2)
    elif what=="mhtSig_m1": calc = lambda jets,jetsUnCorr,corrs,met: calcHTSig(jets,-1)
    elif what=="mhtSig_m2": calc = lambda jets,jetsUnCorr,corrs,met: calcHTSig(jets,-2)
    elif what=="mhtSigST1_p1": calc = lambda jets,jetsUnCorr,corrs,met: calcHTSigST1(jets,jetsUnCorr,met, 1)
    elif what=="mhtSigST1_p2": calc = lambda jets,jetsUnCorr,corrs,met: calcHTSigST1(jets,jetsUnCorr,met, 2)
    elif what=="mhtSigST1_m1": calc = lambda jets,jetsUnCorr,corrs,met: calcHTSigST1(jets,jetsUnCorr,met,-1)
    elif what=="mhtSigST1_m2": calc = lambda jets,jetsUnCorr,corrs,met: calcHTSigST1(jets,jetsUnCorr,met,-2)
    elif what=="mhtSigST2_p1": calc = lambda jets,jetsUnCorr,corrs,met: calcHTSigST2(jets,met,corrs, 1)
    elif what=="mhtSigST2_p2": calc = lambda jets,jetsUnCorr,corrs,met: calcHTSigST2(jets,met,corrs, 2)
    elif what=="mhtSigST2_m1": calc = lambda jets,jetsUnCorr,corrs,met: calcHTSigST2(jets,met,corrs,-1)
    elif what=="mhtSigST2_m2": calc = lambda jets,jetsUnCorr,corrs,met: calcHTSigST2(jets,met,corrs,-2)
    else: 
        print "invalid calc"
        return None

    for i in xrange(tree.GetEntries()):
        tree.GetEntry(i)
        number = getattr(tree, "n"+obj+"Jets")
        rawpt,eta,phi = getattr(tree, obj+"Jets_pt"), getattr(tree, obj+"Jets_eta"), getattr(tree, obj+"Jets_phi")
        met = getattr(tree, obj+post), getattr(tree, obj+post[:-2]+"phi")
        jets = [ ]
        jetsUnCorr = [ ]
        corrs = [] #raw - corrected (to be added to MET)
        for j in xrange(number):
            if abs(eta[j]) > etaCut: continue
            pt = corr.correctedPt(rawpt[j], eta[j]) if corr else rawpt[j]
            if pt > ptCorrCut: 
                jets.append( (pt,eta[j],phi[j]) ) 
                corrs.append( (rawpt[j]-pt,eta[j],phi[j]) ) 
            if rawpt > ptCorrCut:
                jetsUnCorr.append( (rawpt[j],eta[j],phi[j]) )
        ret.append(calc(jets, jetsUnCorr, corrs, met))
    _cache[_key] = ret
    progress.done("done, %d entries" % len(ret))
    return ret
def makeRecoMETArray(tree, what, obj, etaCut):
    if   etaCut <= 1.5: post = "MetBarrel_pt" 
    elif etaCut <= 2.4: post = "MetCentral_pt"
    else:               post = "Met_pt"
    if not tree.GetBranch(obj+post): 
        return None
    progress = _progress("Reading "+obj+post+" ...")
    tree.SetBranchStatus("*",0);
    tree.SetBranchStatus(obj+post,1);
    ret = []
    for i in xrange(tree.GetEntries()):
        tree.GetEntry(i)
        ret.append(getattr(tree,obj+post))
    progress.done("done, %d entries" % len(ret))
    return ret


def genCutCorrArray(corrArray, genArray, genThr):
    if len(genArray) != len(corrArray): raise RuntimeError("Mismatch")
    return [ r for (g,r) in zip(genArray,corrArray) if g > genThr ]

def makeCumulativeHTEff(name, corrArray, xmax, norm=2760.0*11246/1000):
    if len(corrArray) == 0: return None
    nbins = 2000
    ret = ROOT.TH1F("ceff_"+name.replace(" ","_"), "", nbins, 0., xmax);
    for val in corrArray:
        ret.Fill(min(0.9998*xmax, val))
    tot, msum = float(norm)/len(corrArray), 0
    for ib in xrange(0, nbins):
        msum += ret.GetBinContent(nbins-ib) * tot
        ret.SetBinContent(nbins-ib, msum)
    ret.SetDirectory(None)
    return ret

def makeCumulativeHTEffGenCut(name, corrArray, genArray, genThr, xmax, norm):
    return makeCumulativeHTEff(name, genCutCorrArray(corrArray,genArray,genThr), xmax, norm=norm)

def makeROC(effsig,effbkg):
    graph = ROOT.TGraph(effsig.GetNbinsX())
    for i in xrange(1,graph.GetN()+1):
        graph.SetPoint(i-1, effsig.GetBinContent(i), effbkg.GetBinContent(i))
    return graph 

def makeEffHist(name, refArr, corrArr, corrThr, xmax, logxbins=None):
    if logxbins:
        nbins, nratio = int(logxbins[0]), float(logxbins[1])
        if nratio == 1:
            ret = ROOT.TEfficiency(name.replace(" ","_")+"_eff","",nbins,0,xmax)
        else:
            step = pow(nratio, 1.0/nbins)
            base = xmax*(step-1)/(nratio - 1)
            edges = [0]
            for i in xrange(nbins+1):
                edges.append(edges[-1] + base * pow(step, i))
            ret = ROOT.TEfficiency(name.replace(" ","_")+"_eff","",nbins,array('d',edges))
    else:
        ret = ROOT.TEfficiency(name+"_eff","",20,0,xmax)
    for ref,corr in zip(refArr,corrArr):
        ret.Fill( corr >= corrThr , ref )
    ret.SetStatisticOption(ret.kFCP)
    return ret

#from FastPUPPI.NtupleProducer.scripts.respPlots import whats as WHATS
whats = [
    ('l1pfnewpu',[
        ("Puppi",      "L1Puppi$",    ROOT.kGray+1, 20, 1.1,"met",[]), # add var opts
        ("Puppi",      "L1Puppi$",    ROOT.kViolet+1, 20, 1.1,"mht",[1,1]), # add var opts
        ("Puppi",      "L1Puppi$",    ROOT.kOrange-7, 20, 1.1,"mhtST1",[1,1]), # add var opts
        #("Puppi",      "L1Puppi$",    ROOT.kOrange+7, 20, 1.1,"mhtNoJEC",[1,1]), # add var opts
        ("Puppi",      "L1Puppi$",    ROOT.kMagenta-6, 20, 1.1,"mhtST2",[1,1]), # add var opts

        # ("Puppi",      "L1Puppi$",    ROOT.kGreen-3, 20, 1.1,"mhtSig_p1",[1,1]), # add var opts
        # ("Puppi",      "L1Puppi$",    ROOT.kGreen+3, 20, 1.1,"mhtSig_p2",[1,1]), # add var opts
        # ("Puppi",      "L1Puppi$",    ROOT.kRed-4, 20, 1.1,  "mhtSig_m1",[1,1]), # add var opts
        # ("Puppi",      "L1Puppi$",    ROOT.kRed+1, 20, 1.1,  "mhtSig_m2",[1,1]), # add var opts
        # ("Puppi",      "L1Puppi$",    ROOT.kGreen-3, 20, 1.1,"mhtSigST1_p1",[1,1]), # add var opts
        # ("Puppi",      "L1Puppi$",    ROOT.kGreen+3, 20, 1.1,"mhtSigST1_p2",[1,1]), # add var opts
        # ("Puppi",      "L1Puppi$",    ROOT.kRed-4, 20, 1.1,  "mhtSigST1_m1",[1,1]), # add var opts
        # ("Puppi",      "L1Puppi$",    ROOT.kRed+1, 20, 1.1,  "mhtSigST1_m2",[1,1]), # add var opts
        # ("Puppi",      "L1Puppi$",    ROOT.kGreen-3, 20, 1.1,"mhtSigST2_p1",[1,1]), # add var opts
        # ("Puppi",      "L1Puppi$",    ROOT.kGreen+3, 20, 1.1,"mhtSigST2_p2",[1,1]), # add var opts
        # ("Puppi",      "L1Puppi$",    ROOT.kRed-4, 20, 1.1,  "mhtSigST2_m1",[1,1]), # add var opts
        # ("Puppi",      "L1Puppi$",    ROOT.kRed+1, 20, 1.1,  "mhtSigST2_m2",[1,1]), # add var opts

        ("Puppi",      "L1Puppi$",    ROOT.kGreen+3, 20, 1.1,"mhtSig_p2",[1,1]), # add var opts
        ("Puppi",      "L1Puppi$",    ROOT.kGreen-3, 20, 1.1,  "mhtSig_m2",[1,1]), # add var opts
        ("Puppi",      "L1Puppi$",    ROOT.kRed+1, 20, 1.1,"mhtSigST1_p2",[1,1]), # add var opts
        ("Puppi",      "L1Puppi$",    ROOT.kRed-4, 20, 1.1,  "mhtSigST1_m2",[1,1]), # add var opts
        ("Puppi",      "L1Puppi$",    ROOT.kBlue+2, 20, 1.1,"mhtSigST2_p2",[1,1]), # add var opts
        ("Puppi",      "L1Puppi$",    ROOT.kBlue-4, 20, 1.1,  "mhtSigST2_m2",[1,1]), # add var opts
    ]),
]

from optparse import OptionParser
parser = OptionParser("%(prog) infile [ src [ dst ] ]")
parser.add_option("-w", dest="what", default=None, help="Choose set (il1pf, l1pf, ...)")
parser.add_option("-W", dest="what_reg",     default=None, help="Choose set (inputs, l1pf, ...)")
parser.add_option("-P","--plots", dest="plots", default="rate,isorate,roc,effc", help="Choose plot or comma-separated list of plots")
parser.add_option("-j","--jecs", dest="jecs", default="jecs.root", help="Choose JEC file")
parser.add_option("--jm","--jec-method", dest="jecMethod", default="", help="Choose JEC method")
parser.add_option("-R","--raw", dest="rawJets", default=False, action="store_true", help="Don't appy JECs")
parser.add_option("-s", dest="genht",  default=None, type="float", help="Choose gen ht")
parser.add_option("-r", dest="rate",  default="10,20,50", type="string", help="Choose rate [kHz] (for isorate plots, can specify more than one)")
parser.add_option("-l","--label", dest="label",  default=None, type="string", help="Extra label for plots")
parser.add_option("-p", "--pt", dest="pt",  default=30, type="float", help="Choose pt cut")
parser.add_option("-e", "--eta", dest="eta",  default=5.0, type="float", help="Choose eta")
parser.add_option("-v", dest="var",  default="ht", help="Choose variable (ht, met, metCentral, mht, jet<N>, mjj, ptj-mjj<M>)")
parser.add_option("--xlabel","--varlabel", dest="varlabel", default=None, help="X axis label for the variable")
parser.add_option("--xmax", dest="xmax",  default=None, type=float, help="Choose variable")
parser.add_option("--logxbins", dest="logxbins",  default=None, nargs=2, type=float, help="--logxbins N X will make N bins, the last being a factor X larger than the first")
options, args = parser.parse_args()

tfiles = [ROOT.TFile.Open(f) for f in args[:2]]

odir = args[2] 
os.system("mkdir -p "+odir)
os.system("cp %s/src/FastPUPPI/NtupleProducer/python/display/index.php %s/" % (os.environ['CMSSW_BASE'], odir));
ROOT.gROOT.ProcessLine(".x %s/src/FastPUPPI/NtupleProducer/python/display/tdrstyle.cc" % os.environ['CMSSW_BASE']);
ROOT.gStyle.SetOptStat(0)
c1 = ROOT.TCanvas("c1","c1")

ROOT.gSystem.Load("libL1TriggerPhase2L1ParticleFlow")
ROOT.gInterpreter.ProcessLine('#include "L1Trigger/Phase2L1ParticleFlow/src/corrector.h"')

signal, background = [ f.Get("Events") for f in tfiles ]
jecfile = ROOT.TFile.Open(options.jecs)

if options.varlabel is None: options.varlabel = "H_{T}^{miss}"
if options.genht    is None: options.genht    = 150
if options.xmax     is None: options.xmax     = 500
#options.eta = 5.0
qualif = "p_{T}^{corr} > %.0f, |#eta| < %.1f" % (options.pt, options.eta)

for plotkind in options.plots.split(","):
  print "Make plot "+plotkind
  if plotkind != "rate":
    myVar="met"
    genArray = makeGenArray(signal, myVar, options.pt, options.eta)
  rates = map(float, options.rate.split(","))
  if plotkind != "isorate": rates = rates[:1]
  for targetrate in rates:
      for objset,things in whats:
          if options.what and (objset not in options.what.split(",")): continue
          if options.what_reg:
              if not any(re.match(p+"$",objset) for p in options.what_reg.split(",")): 
                  continue
          plots = []
          for name,obj,col,msty,msiz,myVar,myVarOpts in things:
              if "GenAcc$" in obj: continue
              obj = obj.replace("$","")
              if myVar.startswith("met"):
                  jecs = None
              else:
                  jecdirname = obj+"Jets"+( "_"+options.jecMethod if options.jecMethod else "")
                  jecdir = jecfile.GetDirectory(jecdirname)
                  if not jecdir: 
                      print "Missing JECs "+jecdirname+" in "+options.jecs
                      continue
                  jecs = ROOT.l1tpf.corrector(jecdir)
              label = myVar
              if plotkind == "rate":
                  recoArrayB = makeCorrArray(background, myVar, myVarOpts, obj, options.pt, options.eta, jecs)
                  if not recoArrayB: continue
                  plot = makeCumulativeHTEff(name, recoArrayB, options.xmax)
              elif plotkind == "effc":
                  recoArrayS = makeCorrArray(signal, myVar, myVarOpts, obj, options.pt, options.eta, jecs)
                  if not recoArrayS: continue
                  plot = makeCumulativeHTEffGenCut(name, recoArrayS, genArray, options.genht, options.xmax, norm=1)
              elif plotkind == "isorate":
                  recoArrayB = makeCorrArray(background, myVar, myVarOpts, obj, options.pt, options.eta, jecs)
                  if not recoArrayB: continue
                  rateplot = makeCumulativeHTEff(name, recoArrayB, options.xmax)
                  cut = 9999
                  for ix in xrange(1,rateplot.GetNbinsX()+1):
                      if rateplot.GetBinContent(ix) <= targetrate:
                          cut = rateplot.GetXaxis().GetBinLowEdge(ix)
                          break
                  recoArrayS = makeCorrArray(signal, myVar, myVarOpts, obj, options.pt, options.eta, jecs)
                  if not recoArrayS: continue
                  plot = makeEffHist(name, genArray, recoArrayS, cut, options.xmax, logxbins=options.logxbins)
                  label = "%s(%s) > %.0f" % (options.varlabel, myVar,cut)
              elif plotkind == "roc":
                  recoArrayB = makeCorrArray(background, myVar, myVarOpts, obj, options.pt, options.eta, jecs)
                  if not recoArrayB: continue
                  recoArrayS = makeCorrArray(signal,     myVar, myVarOpts, obj, options.pt, options.eta, jecs)
                  if not recoArrayS: continue
                  effsig  = makeCumulativeHTEffGenCut(name+"_s", recoArrayS, genArray, options.genht, options.xmax, norm=1)
                  ratebkg = makeCumulativeHTEff(name+"_b", recoArrayB, options.xmax)
                  if not effsig or not ratebkg: continue
                  plot = makeROC(effsig,ratebkg)
                  msty = 0
              else: raise RuntimeError
              if not plot: continue
              #plot.SetLineWidth(3); plot.SetLineColor(col);  plot.SetMarkerColor(col)
              plot.SetLineWidth(2); plot.SetLineColor(col);  plot.SetMarkerColor(col)
              plot.SetMarkerStyle(msty); plot.SetMarkerSize(msiz)
              plots.append((label,plot))
          if not plots: 
              print "   nothing to plot!"
              continue
          c1.SetLogy(False)
          if plotkind == "rate":
              c1.SetLogy(True)
              frame = ROOT.TH1D("",";L1 %s cut (%s); Minbias rate @ PU200 [kHz]" % (options.varlabel, qualif), 100, 0, options.xmax)
              frame.GetYaxis().SetDecimals(True)
              frame.GetXaxis().SetNdivisions(505)
              frame.GetYaxis().SetRangeUser(0.5, 100e3)
              leg = ROOT.TLegend(0.6,0.99,0.95,0.99-0.06*len(things))
              plotname = '%s%s-%s_eta%s_pt%d' % (options.var, plotkind, objset, options.eta, options.pt)
          elif plotkind == "effc":
              gentext, genpost = "iciency", "" 
              if options.genht > 0:
                  gentext = " (Gen %s > %s)" % (options.varlabel, options.genht)
                  gentpost = "_gen%.0f" % (options.genht)
              frame = ROOT.TH1D("",";L1 %s thresh (%s); Eff%s" % (options.varlabel, qualif, gentext), 100, 0, options.xmax)
              leg = ROOT.TLegend(0.6,0.19,0.95,0.19+0.06*len(things))
              plotname = '%s%s-%s_eta%s_pt%d%s' % (options.var, plotkind, objset, options.eta, options.pt, genpost)
          elif plotkind == "isorate":
              frame = ROOT.TH1D("",";Gen %s (%s); Eff (L1 rate %.0f kHz)" % (options.varlabel, qualif, targetrate), 100, 0, options.xmax)
              frame.GetYaxis().SetDecimals(True)
              leg = ROOT.TLegend(0.65,0.19,0.99,0.19+0.065*len(things))
              plotname = '%s%s-%s_eta%s_pt%d_%.0fkHz' % (options.var, plotkind, objset, options.eta, options.pt, targetrate)
          elif plotkind == "roc":
              c1.SetLogy(True)
              gentext, genpost = "iciency", "" 
              if options.genht > 0:
                  gentext = " (Gen %s > %s)" % (options.varlabel, options.genht)
                  gentpost = "_gen%.0f" % (options.genht)
              frame = ROOT.TH1D("",";Eff%s; Minbias rate @ PU200 [kHz]" % (gentext), 100, 0, 1)
              frame.GetYaxis().SetDecimals(True)
              frame.GetXaxis().SetNdivisions(505)
              frame.GetYaxis().SetRangeUser(0.5, 100e3)
              leg = ROOT.TLegend(0.2,0.98,0.55,0.98-0.06*len(things))
              plotname = '%s%s-%s_eta%s_pt%d%s' % (options.var, plotkind, objset, options.eta, options.pt, genpost)
          frame.Draw()
          if plotkind in ("rate","roc"):
              line = ROOT.TLine()
              line.SetLineStyle(7)
              for y in 40e3, 100, 10:
                  line.DrawLine(frame.GetXaxis().GetXmin(),y,frame.GetXaxis().GetXmax(),y)
          for n,p in plots: 
              p.Draw("PCX SAME" if ("TH1" not in p.ClassName()) else "C SAME")
          for n,p in plots: 
              leg.AddEntry(p, n, "L" if plotkind in ("rate","roc") else "LP")
          leg.Draw()
          
          c1.Print('%s/%s%s.png' % (odir, plotname, ("_"+options.label) if options.label else ""))
          c1.Print('%s/%s%s.pdf' % (odir, plotname, ("_"+options.label) if options.label else ""))
          del frame


