import re, os
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning

from array import array
from math import pow, sqrt
from FastPUPPI.NtupleProducer.scripts.makeJecs import _progress
from FastPUPPI.NtupleProducer.plotTemplate import plotTemplate

ROOT.gROOT.ProcessLine('#include "%s/src/FastPUPPI/NtupleProducer/python/scripts/jetHtSuite.h"' % os.environ['CMSSW_BASE']);

def makeCalcCpp(what, jer):
    if what == "ht": 
        return ROOT.CalcHT()
    if what == "mht": 
        return ROOT.CalcMHT()
    if what == "mhtsig": 
        return ROOT.CalcMHTSig(jer)
    if what.startswith("mhtEvtCorr"):
        return ROOT.CalcMHTCorrEvt(jer, float(what.replace("mhtEvtCorr","")))
    if what.startswith("mhtcorr"):
        return ROOT.CalcMHTCorrJet(jer, float(what.replace("mhtcorr","")))
    if what.startswith("mhttoy"):
        return ROOT.CalcMHTToy(jer, float(what.replace("mhttoy","")))
    if what.startswith("mhtsig-mht"):
        return ROOT.CalcMHTSig_MHTcut(jer, float(what.replace("mhtsig-mht","")))
    if what.startswith("mht-mhtsig"):
        return ROOT.CalcMHT_MHTSigcut(jer, float(what.replace("mht-mhtsig","")))
    if what == "mjj": 
        return ROOT.CalcMJJ()
    if re.match(r"jet\d+$", options.var): 
        return ROOT.CalcJ(int(options.var.replace("jet","")))
    if what.startswith("ptj-mjj"): 
        return ROOT.CalcJ2_MJJcut(float(what.replace("ptj-mjj","")))
    return None

def makeGenArray(tree, what, ptCut, etaCut, jer):
    return makeCorrArray(tree, what, "Gen", ptCut, etaCut, ROOT.nullptr, jer)

def makeCorrArray(tree, what, obj, ptCorrCut, etaCut, corr, jer={}, _cache={}):
    jercalc = jer[obj] if obj in jer else None
    if obj=="Gen" and "L1Puppi" in jer: jercalc = jer["L1Puppi"]
    _key = (id(tree),what,obj,int(ptCorrCut*100),int(etaCut*1000))
    if _key in _cache: return _cache[_key]
    if what == "metmht":
        met = makeCorrArray(tree, "met", obj,     0,      5.0,   corr, jer, _cache=_cache)
        mht = makeCorrArray(tree, "mht", obj, ptCorrCut, etaCut, corr, jer, _cache=_cache)
        ret = ROOT.makeMinimum(met,mht)
        _cache[_key] = ret
        return ret
    if "met" in what:
        ret = makeMETArray(tree, what, obj, etaCut)
        _cache[_key] = ret
        return ret
    if not tree.GetBranch("n"+obj+"Jets"):
        if obj == "Gen": raise RuntimeError("Missing GenJets!");
        return None
    if ('mhtsig' in what) and jercalc==None:
        raise RuntimeError("Missing JER for ",what, obj,"calculation!");
        return None
    cppcalc = makeCalcCpp(what, jercalc)
    if (what.startswith("mhtcorr") or what.startswith("mhttoy") or what.startswith("mhtEvtCorr")) and obj=="Gen":
    # if (what.startswith("mhtcorr") or what.startswith("mhttoy")) and obj=="Gen":
        # for resolution-corrected mht, don't modify gen mht
        cppcalc.SetSigma(0)
    progress = _progress("  Reading "+obj+"Jets in C++...")
    ret = ROOT.makeJetArray(tree, obj, ptCorrCut, etaCut, cppcalc, corr);
    _cache[_key] = ret
    progress.done("done, %d entries" % len(ret))
    return ret

def makeMETArray(tree, what, obj, etaCut):
    if obj == "Gen": obj = "gen" # fix issue with naming convention
    if   etaCut <= 1.5: post = "MetBarrel" 
    elif etaCut <= 2.4: post = "MetCentral"
    else:               post = "Met"
    if not tree.GetBranch(obj+post+"_pt"): 
        if obj == "gen": raise RuntimeError("Missing gen"+post);
        return None
    progress = _progress("  Reading "+obj+post+" in C++...")
    ret = ROOT.makeMetArray(tree, obj+post);
    progress.done("done, %d entries" % len(ret))
    return ret

def makeCumulativeHTEff(name, corrArray, xmax, norm=2760.0*11246/1000):
    return makeCumulativeHTEffGenCut(name, corrArray, None, None, xmax, norm)

def makeCumulativeHTEffGenCut(name, corrArray, genArray, genThr, xmax, norm):
    if len(corrArray) == 0: return None
    nbins = 2000
    ret = ROOT.TH1F("ceff_"+name.replace(" ","_"), "", nbins, 0., xmax);
    if genArray:
        nsel = ROOT.fillTH1FastGenCut(ret, corrArray, genArray, genThr)
    else:
        nsel = ROOT.fillTH1Fast(ret, corrArray)
    if nsel == 0: return None
    tot, msum = float(norm)/nsel, 0
    for ib in xrange(0, nbins):
        msum += ret.GetBinContent(nbins-ib) * tot
        ret.SetBinContent(nbins-ib, msum)
    ret.SetDirectory(ROOT.nullptr)
    return ret

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
    ROOT.fillTEffFast(ret, refArr, corrArr, corrThr)
    ret.SetStatisticOption(ret.kFCP)
    return ret

from FastPUPPI.NtupleProducer.scripts.respPlots import whats as WHATS
whats = WHATS + [
    ('oldcomp',[
        ("Calo",      "L1OldCalo",        ROOT.kViolet+2, 20, 1.5),
        ("TK 5s",     "L1TKV5",           ROOT.kRed+1, 24, 1.5),
        ("PF",        "L1OldPF",          ROOT.kOrange+7, 24, 1.5),
        ("Puppi",     "L1OldPuppi",       ROOT.kBlue+1, 21, 1.5),
        ("Puppi4MET", "L1OldPuppiForMET", ROOT.kAzure+10, 21, 1.5),
    ]),
    ('newcomp',[
        ("Calo",      "L1Calo",        ROOT.kViolet+2, 20, 1.5),
        ("TK 5s",     "L1TKV5",        ROOT.kRed+1, 24, 1.5),
        ("PF",        "L1PF",          ROOT.kOrange+7, 24, 1.5),
        ("Puppi",     "L1Puppi",       ROOT.kBlue+1, 21, 1.5),
        ("Puppi4MET", "L1PuppiForMET", ROOT.kAzure+10, 21, 1.5),
    ]),
    ('l1pfpu_metref',[
        ("Calo",       "L1Calo$",     ROOT.kViolet+1, 21, 1.5),
        ("TK #Deltaz", "L1TKV5$",     ROOT.kGreen+1, 34, 1.2),
        ("RefTK",      "RefL1TrackerEtMiss$",     ROOT.kGreen+3, 34, 1.2),
        ("Puppi",      "L1Puppi$",    ROOT.kRed+1, 20, 1.1),
    ]),
    ('l1pfpu_metnoref',[
        ("Calo",       "L1Calo$",     ROOT.kViolet+1, 21, 1.5),
        ("TK #Deltaz", "L1TKV5$",     ROOT.kGreen+2, 34, 1.2),
        ("Puppi",      "L1Puppi$",    ROOT.kRed+1, 20, 1.1),
    ]),
    ('l1pfpu_metrefonly',[
        ("TK",      "RefL1TrackerEtMiss$",  ROOT.kGreen+2, 34, 1.2),
        ("Puppi",   "L1Puppi$",             ROOT.kRed+1, 20, 1.1),
    ]),
    ('l1pfpu_jetnoref',[
        ("Calo",       "L1Calo$",    ROOT.kViolet+1, 21, 1.5),
        ("TK #Deltaz", "L1TKV5$",    ROOT.kGreen+2, 34, 1.2),
        ("ak4Puppi",   "L1Puppi$",   ROOT.kRed+1, 20, 1.1),
    ]),
    ('l1pfpu_jetref',[
        ("Calo",       "L1Calo$",              ROOT.kViolet+1, 21, 1.5),
        ("RefCalo",    "RefCaloJets$",         ROOT.kViolet+2, 21, 1.5),
        ("TK #Deltaz", "L1TKV5$",              ROOT.kGreen+1, 34, 1.2),
        ("RefTK",      "RefTwoLayerJets$",     ROOT.kGreen+3, 34, 1.2),
        ("ak4Puppi",   "L1Puppi$",             ROOT.kRed+0, 20, 1.1),
        ("RefPuppi",   "RefPhase1PuppiJets$",  ROOT.kRed+2, 20, 0.9), 
    ]),
    ('l1pfpu_jetrefonly',[
        ("Calo",    "RefCaloJets$",         ROOT.kViolet+1, 21, 1.5),
        ("TK",      "RefTwoLayerJets$",     ROOT.kGreen+2, 34, 1.2),
        ("Puppi",   "RefPhase1PuppiJets$",  ROOT.kRed+1, 20, 0.9), 
    ]),
    ('l1pfpu_bitwise',[
        ("CMSSW",       "L1CMSSWPuppi$",      ROOT.kGray+2, 20, 2.0),
        ("Tuned",       "L1Puppi$",           ROOT.kRed+1, 20, 1.7),
        ("Regional",    "L1PuppiRegional$",   ROOT.kViolet+1, 20, 1.4),
        ("BitwisePF",   "L1PuppiBitwisePF$",  ROOT.kAzure+1, 20, 1.1),
        ("BitwisePuppi","L1PuppiBitwise$",    ROOT.kGreen+2, 20, 0.8),
    ]),
]

from optparse import OptionParser
parser = OptionParser("%(prog) infile [ src [ dst ] ]")
parser.add_option("-w", dest="what", default=None, help="Choose set (il1pf, l1pf, ...)")
parser.add_option("-W", dest="what_reg",     default=None, help="Choose set (inputs, l1pf, ...)")
parser.add_option("-P","--plots", dest="plots", default="rate,isorate,roc,effc,plateff,platroc,dist,dist2d", help="Choose plot or comma-separated list of plots") 
parser.add_option("-j","--jecs", dest="jecs", default="jecs.root", help="Choose JEC file")
parser.add_option("--jm","--jec-method", dest="jecMethod", default="", help="Choose JEC method")
parser.add_option("-R","--raw", dest="rawJets", default=False, action="store_true", help="Don't appy JECs")
parser.add_option("-s", dest="genht",  default=None, type="float", help="Choose gen ht")
parser.add_option("-E", dest="eff",  default=None, type="string", help="Choose plateau efficiency")
parser.add_option("-r", dest="rate",  default="10,20,50", type="string", help="Choose rate [kHz] (for isorate plots, can specify more than one)")
parser.add_option("-l","--label", dest="label",  default=None, type="string", help="Extra label for plots")
parser.add_option("-p", "--pt", dest="pt",  default=30, type="float", help="Choose pt cut")
parser.add_option("-e", "--eta", dest="eta",  default=2.4, type="float", help="Choose eta")
parser.add_option("-v", dest="var",  default="ht", help="Choose variable (ht, met, metCentral, mht, jet<N>, mjj, ptj-mjj<M>, mhtsig, mhtsig-mht<M>, mht-mhtsig<M>, mhtcorr<M>, mhttoy<M>)")
parser.add_option("--xlabel","--varlabel", dest="varlabel", default=None, help="X axis label for the variable")
parser.add_option("--xmax", dest="xmax",  default=None, type=float, help="Choose variable")
parser.add_option("--logxbins", dest="logxbins",  default=None, nargs=2, type=float, help="--logxbins N X will make N bins, the last being a factor X larger than the first")
parser.add_option("--print", dest="printQualityPlots", action="store_true", default=False, help="Make print-quality plots (incl. pdf & eps)")
parser.add_option("--jer", dest="jer", default="", help="Choose jet energy resolution (JER) file")
options, args = parser.parse_args()

tfiles = [ROOT.TFile.Open(f) for f in args[:2]]

odir = args[2] 
plotter = plotTemplate(odir, defaultExts = (["png","eps","pdf"] if options.printQualityPlots else ["png"]))

ROOT.gSystem.Load("libL1TriggerPhase2L1ParticleFlow")
ROOT.gInterpreter.ProcessLine('#include "L1Trigger/Phase2L1ParticleFlow/src/corrector.h"')

if options.var == "ht":
    if options.varlabel is None: options.varlabel = "H_{T}"
    if options.genht    is None: options.genht    = 300
    if options.xmax     is None: options.xmax     = 1000
    if options.eff      is None: options.eff      = "0.5,0.9,0.95"
    qualif = "p_{T}^{corr} > %.0f, |#eta| < %.1f" % (options.pt, options.eta)
    what = options.var
elif options.var == "mht":
    if options.varlabel is None: options.varlabel = "H_{T}^{miss}"
    if options.genht    is None: options.genht    = 150
    if options.xmax     is None: options.xmax     = 500
    if options.eff      is None: options.eff      = "0.5,0.9,0.95"
    qualif = "p_{T}^{corr} > %.0f, |#eta| < %.1f" % (options.pt, options.eta)
    what = options.var
elif options.var.startswith("mhtcorr") or options.var.startswith("mhtEvtCorr"):
    if options.varlabel is None: options.varlabel = "H_{T}^{miss}(1-x#sigma)"
    if options.genht    is None: options.genht    = 150
    if options.xmax     is None: options.xmax     = 500
    if options.eff      is None: options.eff      = "0.95"
    qualif = "p_{T}^{corr} > %.0f, |#eta| < %.1f" % (options.pt, options.eta)
    what = options.var
elif options.var.startswith("mhttoy"):
    if options.varlabel is None: options.varlabel = "H_{T}^{miss}(toy, X pct)"
    if options.genht    is None: options.genht    = 150
    if options.xmax     is None: options.xmax     = 500
    if options.eff      is None: options.eff      = "0.95"
    qualif = "p_{T}^{corr} > %.0f, |#eta| < %.1f" % (options.pt, options.eta)
    what = options.var
elif options.var == "mhtsig":
    if options.varlabel is None: options.varlabel = "H_{T}^{miss} significance"
    if options.genht    is None: options.genht    = 150
    if options.xmax     is None: options.xmax     = 15
    if options.eff      is None: options.eff      = "0.5,0.9,0.95"
    qualif = "p_{T}^{corr} > %.0f, |#eta| < %.1f" % (options.pt, options.eta)
    what = options.var
elif options.var.startswith("mhtsig-mht"):
    if options.varlabel is None: options.varlabel = "H_{T}^{miss} significance"
    if options.genht    is None: options.genht    = 150
    if options.xmax     is None: options.xmax     = 15
    if options.eff      is None: options.eff      = "0.5,0.9,0.95"
    qualif = "p_{T}^{corr} > %.0f, |#eta| < %.1f, mH_{T}> %s" % (options.pt, options.eta, options.var.replace("mhtsig-mht",""))
    what = options.var
elif options.var.startswith("mht-mhtsig"):
    if options.varlabel is None: options.varlabel = "H_{T}^{miss}"
    if options.genht    is None: options.genht    = 150
    if options.xmax     is None: options.xmax     = 500
    if options.eff      is None: options.eff      = "0.5,0.9,0.95"
    qualif = "p_{T}^{corr} > %.0f, |#eta| < %.1f, #sigma(mH_{T})> %s" % (options.pt, options.eta, options.var.replace("mht-mhtsig",""))
    what = options.var
elif options.var == "metmht":
    if options.varlabel is None: options.varlabel = "E_{T}^{miss} - H_{T}^{miss}"
    if options.genht    is None: options.genht    = 150
    if options.xmax     is None: options.xmax     = 500
    if options.eff      is None: options.eff      = "0.5,0.9,0.95"
    qualif = "p_{T}^{corr} > %.0f, |#eta| < %.1f" % (options.pt, options.eta)
    what = options.var
elif options.var == "mjj":
    if options.varlabel is None: options.varlabel = "m(jj)"
    if options.genht    is None: options.genht    = 750
    if options.xmax     is None: options.xmax     = 2000
    if options.eff      is None: options.eff      = "0.9,0.95"
    qualif = "p_{T}^{corr} > %.0f, |#eta| < %.1f" % (options.pt, options.eta)
    what = options.var
elif re.match(r"jet\d+$", options.var):
    ijet = int(options.var.replace("jet",""))
    if options.varlabel is None: options.varlabel = "jet %d p_{T}" % ijet
    if options.genht    is None: options.genht    = 200 if ijet <= 2 else ( 45 if ijet <= 4 else  30)
    if options.xmax     is None: options.xmax     = 400 if ijet <= 2 else (150 if ijet <= 4 else 100)
    if options.eff      is None: options.eff      = "0.95" if ijet <= 2 else "0.9,0.95"
    options.pt = 10
    qualif = "|#eta| < %.1f" % (options.eta)
    what = options.var
elif options.var.startswith("ptj-mjj"):
    if options.varlabel is None: options.varlabel = "jet p_{T}"
    if options.genht    is None: options.genht    = 100
    if options.xmax     is None: options.xmax     = 300
    if options.eff      is None: options.eff      = "0.9,0.95"
    options.pt = 20 # slightly increase the cut, to avoid too large combinatoric.
    qualif = "|#eta| < %.1f, m(jj) > %s" % (options.eta, options.var.replace("ptj-mjj",""))
    what = options.var
elif options.var.startswith("met"):
    if options.varlabel is None: options.varlabel = "E_{T}^{miss}"
    if options.genht    is None: options.genht    = 150
    if options.xmax     is None: options.xmax     = 500
    if options.eff      is None: options.eff      = "0.5,0.9,0.95"
    what = "met"
    options.eta = 5.0 
    if "Central" in options.var:  options.eta = 2.4 
    elif "Barrel" in options.var: options.eta = 1.5
    qualif = "|#eta| < %.1f" % options.eta
else:
    raise RuntimeError("Unknown variable "+options.var)

signal, background = [ f.Get("Events") for f in tfiles ]
jecfile = ROOT.TFile.Open(options.jecs)

if ("mhtsig" in options.var) and len(options.jer)==0:
    exit("Significance calc requires --jer option")
jerTool={}
if options.jer:
    print "Loading JER Tool"
    jerfile = ROOT.TFile.Open(options.jer)
    fits={}
    for k in jerfile.GetListOfKeys():
        kname = k.GetName()
        if '_jet_eta_' in kname:
            kind, etabin = kname.split('_jet_eta_')
            if not(kind in fits): fits[kind]=[]
            fits[kind].append( (etabin, jerfile.Get(kname)) )
            elow, ehi = map(lambda x:float(x)/10., etabin.split('_'))
    for kind in fits:
        fits[kind].sort()
        bins = [x[0].split('_')[0] for x in fits[kind]] + [ fits[kind][-1][0].split('_')[1] ]
        bins = array('d', map(lambda x : float(x)/10., bins))
        funcs = ROOT.vector('TF1')(0)
        for x in fits[kind]: funcs.push_back(x[1])
        # print kind, bins, funcs
        jerTool[kind] = ROOT.JetResolutionCalc(len(bins)-1, bins, funcs)
        print "For ",kind," cfg, loaded JER abs(eta) bins:", bins.tolist()
    

def makePlatEffPlot(signal, background, what, obj, ptcut, jecs, jer, plotparam, _cache={}):
    _key = (id(signal),id(background),what,obj,str(ptcut),str(plotparam))
    if _key in _cache: 
        #print "  retrieved plateff for %s, %s, %s from cache" % (what, obj, plotparam)
        return _cache[_key]
    # ok there we go
    recoArrayB = makeCorrArray(background, what, obj, ptcut, options.eta, jecs, jer)
    if not recoArrayB: return (None, None)
    rateplot = makeCumulativeHTEff(name, recoArrayB, options.xmax)
    recoArrayS = makeCorrArray(signal, what, obj, ptcut, options.eta, jecs, jer)
    if not recoArrayS: return (None, None)
    #platprogress = _progress("  making plateff for %s, %s, %s" % (what, obj, plotparam))
    def effForRate(rate):
      cut = 9999
      for ix in xrange(1,rateplot.GetNbinsX()+1):
          if rateplot.GetBinContent(ix) <= rate:
              cut = rateplot.GetXaxis().GetBinLowEdge(ix)
              break
      #print "Cut for %s @ rate %g: %g" % (what+obj, rate, cut)
      plot = makeEffHist(name, genArray, recoArrayS, cut, options.xmax, logxbins=options.logxbins)
      eff = plot.GetEfficiency(plot.FindFixBin(options.genht))
      #print "Eff at %g: %g" % (options.genht, eff)
      return eff, plot, cut
    rate = 50
    while True:
      eff, plot, cut = effForRate(rate)
      if not eff: break
      if eff < plotparam:
          if rate >= 30e3: 
              eff = None; break
          rate = min(rate*5, 30e3)
      else:
          break
    if not eff: return (None, None)
    #print "Upper bound: eff %g at rate %g" % (eff, rate) 
    maxrate = rate; rate = rate / 5;
    while True:
      eff, plot, cut = effForRate(rate)
      if not eff: break
      if eff > plotparam:
          rate /= 2; 
      else:
          break
    if not eff: return (None, None)                 
    # print "Lower bound: eff %g at rate %g" % (eff, rate) 
    minrate = rate
    while True:
      rate = sqrt(maxrate * minrate)
      eff, plot, cut = effForRate(rate)
      if not eff: 
          # print "Bisection failed?"
          break
      if eff < plotparam:
          # print "New lower bound: eff %g at rate %g" % (eff, rate)
          minrate = rate;
      else:
          # print "New upper bound: eff %g at rate %g" % (eff, rate)
          maxrate = rate
      if maxrate/minrate < 1.1: 
          break
    #platprogress.done()
    if not eff: return (None, None)
    if rate < 10:    ratestr = "%.1fkHz" % rate
    elif rate < 500: ratestr = "%.0fkHz" % rate
    else:            ratestr = "%.1fMHz" % (rate/1000)
    label = "%s @ %s" % (name, ratestr)
    _cache[_key] = (plot,label)
    return (plot,label)

def makePlatRocPlot(signal, background, what, obj, ptcut, jecs, jer, plotparam, _cache={}):
    _key = (id(signal),id(background),what,obj,str(ptcut),str(plotparam))
    if _key in _cache: 
        #print "  retrieved platroc for %s, %s, %s from cache" % (what, obj, plotparam)
        return _cache[_key]
    # ok there we go
    #platprogress = _progress("  making platroc for %s, %s, %s" % (what, obj, plotparam))
    recoArrayB = makeCorrArray(background, what, obj, ptcut, options.eta, jecs, jer)
    if not recoArrayB: return (None,None)
    recoArrayS = makeCorrArray(signal, what, obj, ptcut, options.eta, jecs, jer)
    if not recoArrayS: return (None,None)
    rateplot = makeCumulativeHTEff(name, recoArrayB, options.xmax)
    def platForRate(rate):
      cut = None
      for ix in xrange(1,rateplot.GetNbinsX()+1):
          if rateplot.GetBinContent(ix) <= rate:
              cut = rateplot.GetXaxis().GetBinLowEdge(ix)
              break
      if cut is None:
          return None
      #print "Cut for %s @ rate %g: %g" % (what+obj, rate, cut)
      plot = makeEffHist(name, genArray, recoArrayS, cut, options.xmax, logxbins=options.logxbins)
      hist = plot.GetTotalHistogram(); xaxis = hist.GetXaxis()
      for i in xrange(2,hist.GetNbinsX()+1):
          if plot.GetEfficiency(i) > plotparam:
              graph = plot.CreateGraph()
              xmin, xmax = xaxis.GetBinLowEdge(i-1), xaxis.GetBinCenter(i)
              if graph.Eval(xmin,ROOT.nullptr,"S") > plotparam: 
                  #print "ERROR xmin for %s @ %g (%g): i = %d" % (what+obj, rate, cut, i)
                  pass
              elif graph.Eval(xmax,ROOT.nullptr,"S") < plotparam:
                  #print "ERROR xmax for %s @ %g (%g): i = %d" % (what+obj, rate, cut, i)
                  pass
              else:
                  xmid = 0.5*(xmax+xmin)
                  while abs(xmax-xmin) > 0.05*(abs(xmax)+abs(xmin)):
                      if graph.Eval(xmid,ROOT.nullptr,"S") < plotparam:
                          xmin = xmid
                      else:
                          xmax = xmid
                      xmid = 0.5*(xmax+xmin)
                  #print "Plateau at %g eff for %s @ %g: %g" % (plotparam, what+obj, rate, xmid)
                  return xmid
              break
      return None
    points = []
    rate = 1
    while rate <= 30e3:
      plat = platForRate(rate)
      if plat: 
          if len(points) == 0 or points[-1][0] != plat:
              points.append((plat,rate))
      rate *= 1.2
    #platprogress.done(" done, with %d points" % len(points))
    if not points: return (None,None)
    plot = ROOT.TGraph(len(points))
    for i,(x,y) in enumerate(points):
      plot.SetPoint(i,x,y)
    label = name
    _cache[_key] = (plot,label)
    return (plot,label)

print "Plotting for %s (%s)" % (options.var, options.varlabel)
for plotkind in options.plots.split(","):
  progress = _progress("Make plot %s: \n" % plotkind)
  if plotkind != "rate":
    genArray = makeGenArray(signal, what, options.pt, options.eta, jerTool)
  if plotkind == "isorate":
      plotparams = map(float, options.rate.split(","))
  elif plotkind in ("platroc", "plateff"):
      plotparams = map(float, options.eff.split(","))
  elif plotkind == "dist":
      plotparams = ["signal","signalGen","bkg"]
  elif plotkind == "dist2d":
      plotparams = ["signal","bkg"]
  else:
      plotparams = [None]
  for plotparam in plotparams:
      for objset,things in whats:
          if (plotkind == "dist2d") and not (("mhtsig" in options.var) or ("mhtcorr" in options.var) or ("mhttoy" in options.var)): continue
          if options.what and (objset not in options.what.split(",")): continue
          if options.what_reg:
              if not any(re.match(p+"$",objset) for p in options.what_reg.split(",")): 
                  continue
          plots = []
          for name,obj,col,msty,msiz in things:
              if "GenAcc$" in obj: continue
              obj = obj.replace("$","")
              if options.var.startswith("met") or obj.startswith("Ref"):
                  jecs = None
              else:
                  jecdirname = obj+"Jets"+( "_"+options.jecMethod if options.jecMethod else "")
                  jecdir = jecfile.GetDirectory(jecdirname)
                  if not jecdir: 
                      print "Missing JECs "+jecdirname+" in "+options.jecs
                      continue
                  jecs = ROOT.l1tpf.corrector(jecdir)
              label = name
              ptcut = options.pt
              if "RefTwoLayerJets" in obj: ptcut = 5
              if plotkind == "rate":
                  recoArrayB = makeCorrArray(background, what, obj, ptcut, options.eta, jecs, jerTool)
                  if not recoArrayB: continue
                  plot = makeCumulativeHTEff(name, recoArrayB, options.xmax)
              elif plotkind == "effc":
                  recoArrayS = makeCorrArray(signal, what, obj, ptcut, options.eta, jecs, jerTool)
                  if not recoArrayS: continue
                  plot = makeCumulativeHTEffGenCut(name, recoArrayS, genArray, options.genht, options.xmax, norm=1)
              elif plotkind == "isorate":
                  targetrate = plotparam
                  recoArrayB = makeCorrArray(background, what, obj, ptcut, options.eta, jecs, jerTool)
                  if not recoArrayB: continue
                  rateplot = makeCumulativeHTEff(name, recoArrayB, options.xmax)
                  if not rateplot: continue
                  cut = 9999
                  for ix in xrange(1,rateplot.GetNbinsX()+1):
                      if rateplot.GetBinContent(ix) <= targetrate:
                          cut = rateplot.GetXaxis().GetBinLowEdge(ix)
                          break
                  recoArrayS = makeCorrArray(signal, what, obj, ptcut, options.eta, jecs, jerTool)
                  if not recoArrayS: continue
                  plot = makeEffHist(name, genArray, recoArrayS, cut, options.xmax, logxbins=options.logxbins)
                  label = "%s(%s) > %.0f" % (options.varlabel, name,cut)
              elif plotkind == "roc":
                  recoArrayB = makeCorrArray(background, what, obj, ptcut, options.eta, jecs, jerTool)
                  if not recoArrayB: continue
                  recoArrayS = makeCorrArray(signal,     what, obj, ptcut, options.eta, jecs, jerTool)
                  if not recoArrayS: continue
                  effsig  = makeCumulativeHTEffGenCut(name+"_s", recoArrayS, genArray, options.genht, options.xmax, norm=1)
                  ratebkg = makeCumulativeHTEff(name+"_b", recoArrayB, options.xmax)
                  if not effsig or not ratebkg: continue
                  plot = ROOT.makeROCFast(effsig,ratebkg)
                  msty = 0
              elif plotkind == "plateff":
                  (plot, label) = makePlatEffPlot(signal, background, what, obj, ptcut, jecs, jerTool, plotparam)
              elif plotkind == "platroc":
                  (plot, label) = makePlatRocPlot(signal, background, what, obj, ptcut, jecs, jerTool, plotparam)
              elif 'dist' in plotkind:
                  if plotparam=="signal":
                      arr = makeCorrArray(signal,     what, obj, ptcut, options.eta, jecs, jerTool)
                  elif plotparam=="bkg":
                      arr = makeCorrArray(background, what, obj, ptcut, options.eta, jecs, jerTool)
                  elif plotparam=="signalGen":
                      arr = genArray
                  label = plotkind+"_"+plotparam+"_"+name.replace(" ","_")
                  if plotkind == "dist":
                      print '1d', label
                      plot = ROOT.TH1F(label, "", 100, 0., options.xmax);
                      nsel = ROOT.fillTH1Fast(plot, arr)
                  elif plotkind == "dist2d":
                      if plotparam=="signal":
                          arrMHT = makeCorrArray(signal,     "mht", obj, ptcut, options.eta, jecs, jerTool)
                      elif plotparam=="bkg":
                          arrMHT = makeCorrArray(background, "mht", obj, ptcut, options.eta, jecs, jerTool)
                      print '2d', label
                      plot = ROOT.TH2F(label, "", 80, 0., 300, 80, 0., options.xmax);
                      nsel = ROOT.fillTH2Fast(plot, arrMHT, arr)
                  if plot.GetMaximum(): plot.Scale(1./plot.GetMaximum())
              else: raise RuntimeError
              if not plot: continue
              plot.SetLineWidth(3); plot.SetLineColor(col);  plot.SetMarkerColor(col)
              plot.SetMarkerStyle(msty); plot.SetMarkerSize(msiz)
              plots.append((label,plot))
          if not plots: 
              print "   nothing to plot!"
              continue
          plotter.SetLogy(False)
          if plotkind == "rate":
              plotter.SetLogy(True)
              frame = ROOT.TH1D("",";L1 %s cut (%s); Minbias rate @ PU200 [kHz]" % (options.varlabel, qualif), 100, 0, options.xmax)
              frame.GetYaxis().SetDecimals(True)
              frame.GetXaxis().SetNdivisions(505)
              frame.GetYaxis().SetRangeUser(0.5, 100e3)
              leg = ROOT.TLegend(0.56,0.93,0.93,0.93-0.055*len(things))
              plotname = '%s%s-%s_eta%s_pt%d' % (options.var, plotkind, objset, options.eta, options.pt)
          elif plotkind == "effc":
              gentext, genpost = "iciency", "" 
              if options.genht > 0:
                  gentext = " (Gen %s > %s)" % (options.varlabel, options.genht)
                  genpost = "_gen%.0f" % (options.genht)
              frame = ROOT.TH1D("",";L1 %s thresh (%s); Eff%s" % (options.varlabel, qualif, gentext), 100, 0, options.xmax)
              leg = ROOT.TLegend(0.6,0.93,0.93,0.93-0.055*len(things))
              plotname = '%s%s-%s_eta%s_pt%d%s' % (options.var, plotkind, objset, options.eta, options.pt, genpost)
          elif plotkind == "isorate":
              frame = ROOT.TH1D("",";Gen %s (%s); Eff (L1 rate %.0f kHz)" % (options.varlabel, qualif, plotparam), 100, 0, options.xmax)
              frame.GetYaxis().SetDecimals(True)
              leg = ROOT.TLegend(0.50,0.18,0.93,0.18+0.055*len(things))
              plotname = '%s%s-%s_eta%s_pt%d_%.0fkHz' % (options.var, plotkind, objset, options.eta, options.pt, plotparam)
          elif plotkind == "plateff":
              frame = ROOT.TH1D("",";Gen %s (%s); Eff" % (options.varlabel, qualif), 100, 0, options.xmax)
              frame.GetYaxis().SetDecimals(True)
              leg = ROOT.TLegend(0.30,0.18,0.93,0.18+0.055*len(things))
              plotname = '%s%s-%s_eta%s_pt%d_eff%.2fat%g' % (options.var, plotkind, objset, options.eta, options.pt, plotparam, options.genht)
          elif plotkind == "roc":
              plotter.SetLogy(True)
              gentext, genpost = "iciency", "" 
              if options.genht > 0:
                  gentext = " (Gen %s > %s)" % (options.varlabel, options.genht)
                  genpost = "_gen%.0f" % (options.genht)
              frame = ROOT.TH1D("",";Eff%s; Minbias rate @ PU200 [kHz]" % (gentext), 100, 0, 1)
              frame.GetYaxis().SetDecimals(True)
              frame.GetXaxis().SetNdivisions(505)
              frame.GetYaxis().SetRangeUser(0.5, 100e3)
              leg = ROOT.TLegend(0.2,0.93,0.55,0.93-0.055*len(things))
              plotname = '%s%s-%s_eta%s_pt%d%s' % (options.var, plotkind, objset, options.eta, options.pt, genpost)
          elif plotkind == "platroc":
              plotter.SetLogy(True)
              xtit = "Gen threshold at %g%% eff" % (plotparam*100)
              frame = ROOT.TH1D("",";%s; Minbias rate @ PU200 [kHz]" % (xtit), 100, 0, options.xmax)
              frame.GetYaxis().SetDecimals(True)
              #frame.GetXaxis().SetNdivisions(505)
              frame.GetYaxis().SetRangeUser(0.5, 100e3)
              leg = ROOT.TLegend(0.6,0.93,0.93,0.93-0.055*len(things))
              plotname = '%s%s-%s_eta%s_pt%d_eff%.3f' % (options.var, plotkind, objset, options.eta, options.pt, plotparam)
          elif plotkind == "dist":
              plotter.SetLogy(False)
              frame = ROOT.TH1D("",";%s %s; arbitrary" % (plotparam, options.varlabel), 100, 0, options.xmax)
              frame.GetYaxis().SetDecimals(True)
              frame.GetXaxis().SetNdivisions(505)
              frame.GetYaxis().SetRangeUser(0., 1.2)
              leg = ROOT.TLegend(0.2,0.93,0.55,0.93-0.055*len(things))
              plotname = '%s%s-%s_eta%s_pt%d%s' % (options.var, plotkind, objset, options.eta, options.pt, plotparam)
          elif plotkind == "dist2d":
              plotter.SetLogy(False)
              plotter.canvas.SetLogz(True)
              plotter.canvas.SetRightMargin(0.15)
              frame = ROOT.TH2D("",";MHT [GeV];%s %s; arbitrary" % (plotparam, options.varlabel), 80,0,300, 80, 0, options.xmax)
              frame.GetYaxis().SetDecimals(True)
              frame.GetXaxis().SetNdivisions(505)
              #frame.GetZaxis().SetRangeUser(0., 1.2)
              frame.GetZaxis().SetRangeUser(1e-6, 1.5)
              leg = ROOT.TLegend(0.2,0.93,0.55,0.93-0.055*len(things))
              plotname = '%s%s-%s_eta%s_pt%d%s' % (options.var, plotkind, objset, options.eta, options.pt, plotparam)
          frame.Draw()
          if plotkind in ("rate","roc","platroc"):
              line = ROOT.TLine()
              line.SetLineStyle(7)
              for y in 40e3, 100, 10:
                  line.DrawLine(frame.GetXaxis().GetXmin(),y,frame.GetXaxis().GetXmax(),y)
          for n,p in plots: 
              if plotkind == "platroc":
                  p.Draw("PLX SAME")
              elif plotkind == "dist2d":
                  p.Draw("COLZ SAME")
              else:
                  p.Draw("PCX SAME" if ("TH1" not in p.ClassName()) else "C SAME")
          for n,p in plots: 
              leg.AddEntry(p, n, "L" if plotkind in ("rate","roc") else "LP")
          if plotkind != "dist2d": leg.Draw()
          plotter.decorations()
          #printprogress = _progress("  Printing plot %s%s (plot)" % (plotname, ("_"+options.label) if options.label else ""))
          plotter.Print('%s%s' % (plotname, ("_"+options.label) if options.label else ""))
          fout = ROOT.TFile.Open('%s/%s%s.root' % (odir, plotname, ("_"+options.label) if options.label else ""), "RECREATE")
          fout.WriteTObject(frame,"frame")
          for n,p in plots: 
              p.SetTitle(n)
              fout.WriteTObject(p)
          fout.Close()
          del frame
          #printprogress.done()
  progress.done("  Done plot %s" % plotkind)
  print ""


