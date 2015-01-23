import math
import ROOT

ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)
#ROOT.gStyle.SetFillStyle(ROOT.kWhite)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetFrameBorderMode(ROOT.kWhite)
ROOT.gStyle.SetFrameFillColor(ROOT.kWhite)
ROOT.gStyle.SetCanvasBorderMode(ROOT.kWhite)
ROOT.gStyle.SetCanvasColor(ROOT.kWhite)
ROOT.gStyle.SetPadBorderMode(ROOT.kWhite)
ROOT.gStyle.SetPadColor(ROOT.kWhite)
ROOT.gStyle.SetStatColor(ROOT.kWhite)
ROOT.gStyle.SetErrorX(0)

canvas = ROOT.TCanvas('canvas', '', 100, 100, 1200, 800)
canvas.SetGridx()
canvas.SetGridy()

file_out = ROOT.TFile('histograms.root','RECREATE')

collection_suffixes = []
collection_suffixes.append('ID'       )
collection_suffixes.append('isolation')
collection_suffixes.append('total'    )

cut_suffixes = {}
for cs in collection_suffixes:
    cut_suffixes[cs] = []

class var_object:
    def __init__(self, name, nBins, lower, upper, xaxis, unit):
        self.name  = name
        self.nBins = nBins
        self.lower = lower
        self.upper = upper
        self.xaxis = xaxis
        self.unit  = unit
        self.per   = 1.0*(self.upper-self.lower)/self.nBins
        self.histogram = ROOT.TH1F('hBase_%s'%self.name,'',self.nBins,self.lower,self.upper)
        xaxis_title = self.xaxis if self.unit=='' else '%s [%s]'%(self.xaxis,self.unit)
        self.histogram.GetXaxis().SetTitle(xaxis_title)
        self.histogram.GetYaxis().SetTitle('electrons per %.4f %s'%(self.per,self.unit))
        self.histogram.SetFillColor(ROOT.kRed)
        self.histogram.SetLineColor(ROOT.kBlack)
        self.histogram.SetMarkerColor(ROOT.kBlack)
        self.histogram.SetMarkerStyle(20)

var_names = []
var_names.append('Et'             )
var_names.append('eta'            )
var_names.append('EcalDriven'     )
var_names.append('dEtaIn'         )
var_names.append('dPhiIn'         )
var_names.append('HOverE'         )
var_names.append('SigmaIetaIeta'  )
var_names.append('E1x5OverE5x5'   )
var_names.append('E2x5OverE5x5'   )
var_names.append('missingHits'    )
var_names.append('dxyFirstPV'     )
var_names.append('isolEMHadDepth1')
var_names.append('IsolPtTrks'     )

vars = {}
vars['Et'             ] = var_object('Et'             ,  100,     0, 3000, 'E_{T}(e)'                , 'GeV')
vars['eta'            ] = var_object('eta'            ,  100,    -3,    3, '#eta_{T}(e)'             , ''   )
vars['EcalDriven'     ] = var_object('EcalDriven'     ,    3,  -1.5,  1.5, 'isEcalDriven(e)'         , ''   )
vars['dEtaIn'         ] = var_object('dEtaIn'         ,  100, -0.05, 0.05, '#Delta#eta_{in}(e)'      , ''   )
vars['dPhiIn'         ] = var_object('dPhiIn'         ,  100,  -0.1,  0.1, '#Delta#phi_{in}(e)'      , ''   )
vars['HOverE'         ] = var_object('HOverE'         ,  100,     0,    1, 'H/E(e)'                  , ''   )
vars['SigmaIetaIeta'  ] = var_object('SigmaIetaIeta'  ,  100,     0, 0.05, '#sigma_{i#etai#eta}(e)'  , ''   )
vars['E1x5OverE5x5'   ] = var_object('E1x5OverE5x5'   ,  100,     0,    1, 'E_{1x5}/E_{5x5}(e)'      , ''   )
vars['E2x5OverE5x5'   ] = var_object('E2x5OverE5x5'   ,  100,   0.5,    1, 'E_{2x5}/E_{5x5}(e)'      , ''   )
vars['missingHits'    ] = var_object('missingHits'    ,   11,  -0.5, 10.5, 'missing hits (e)'        , ''   )
vars['dxyFirstPV'     ] = var_object('dxyFirstPV'     ,  100,  -0.1,  0.1, 'd_{xy}^{PV}(e)'          , 'mm' )
vars['isolEMHadDepth1'] = var_object('isolEMHadDepth1',  100,     0,   50, 'isolation_{HadDepth1}(e)', 'GeV')
vars['IsolPtTrks'     ] = var_object('IsolPtTrks'     ,  100,     0,   10, 'isolation_{Tracks}(e)'   , ''   )

cut_suffixes['ID'].append('Et'           )
cut_suffixes['ID'].append('eta'          )
cut_suffixes['ID'].append('EcalDriven'   )
cut_suffixes['ID'].append('dEtaIn'       )
cut_suffixes['ID'].append('dPhiIn'       )
cut_suffixes['ID'].append('HOverE'       )
cut_suffixes['ID'].append('SigmaIetaIeta')
cut_suffixes['ID'].append('E1x5OverE5x5' )
cut_suffixes['ID'].append('E2x5OverE5x5' )
cut_suffixes['ID'].append('missingHits'  )
cut_suffixes['ID'].append('dxyFirstPV'   )

cut_suffixes['isolation'].append('isolEMHadDepth1')
cut_suffixes['isolation'].append('IsolPtTrks'     )

cut_suffixes['total'] = cut_suffixes['ID']+cut_suffixes['isolation']

cutflow_prefixes = []
cutflow_prefixes.append('HEEP_cutflow41'     )
cutflow_prefixes.append('HEEP_cutflow50_25ns')
cutflow_prefixes.append('HEEP_cutflow50_50ns')

regions = ['barrel','endcap']

class cutflow_instance_object:
    def __init__(self, name, value, success):
        self.name    = name
        self.value   = value
        self.success = success

class electron_object:
    def __init__(self, cutflows_in, tree, i_el):
        self.isValid = i_el < len(tree.HEEP_cutflow41_total)
        if self.isValid == False:
            return
        self.pt     = tree.gsf_pt[i_el]
        self.eta    = tree.gsf_eta[i_el]
        self.phi    = tree.gsf_phi[i_el]
        self.energy = tree.gsf_energy[i_el]
        self.region = 'barrel' if abs(self.eta) < 1.4442 else 'endcap'
        self.p4 = ROOT.TLorentzVector()
        self.p4.SetPtEtaPhiE(self.pt, self.eta, self.phi, self.energy)
        
        self.DR_matched = False
        self.DR_best_energy_ratio = 0
        for i_mc in range(0,len(tree.mc_energy)):
            mcp4 = ROOT.TLorentzVector()
            mcp4.SetPtEtaPhiE( tree.mc_pt[i_mc], tree.mc_eta[i_mc], tree.mc_phi[i_mc], tree.mc_energy[i_mc] )
            if self.p4.DeltaR(mcp4) < 0.15:
                self.DR_matched = True
                if self.p4.E()/mcp4.E() > self.DR_best_energy_ratio:
                    self.DR_best_energy_ratio = self.p4.E()/mcp4.E()
                    self.mcp4 = mcp4
        if self.DR_best_energy_ratio < 0.9:
            self.isValid = False
        
        self.cutflows = {}
        self.cutflow_bitmaps = {}
        for c in cutflows_in:
            cname = c.name
            self.cutflows[cname] = []
            self.cutflow_bitmaps[cname] = 1
            self.cutflows[c.name].append(cutflow_instance_object('raw',-999,True))
            for i_cn in range(0,len(c.cut_names)):
                cn = c.cut_names[i_cn]
                value   = getattr(tree, '%s_value'%cn)[i_el]
                success = getattr(tree, cn           )[i_el]
                self.cutflows[cname].append(cutflow_instance_object(cn,value,success))
                if success:
                    self.cutflow_bitmaps[cname] += math.pow(2,i_cn+1)
            self.cutflow_bitmaps[cname] = int(self.cutflow_bitmaps[cname])

class cutflow:
    def __init__(self, prefix, suffix):
        self.prefix = prefix
        self.suffix = suffix
        self.name   = '%s_%s'%(self.prefix,self.suffix)
        self.cut_names = []
        self.nEvents               = {}
        self.NM1_nEvents           = {}
        self.cumulative_nEvents    = {}
        self.nElectrons            = {}
        self.NM1_nElectrons        = {}
        self.cumulative_nElectrons = {}
        self.cut_suffixes = cut_suffixes[self.suffix]
        for cs in cut_suffixes[self.suffix]:
            self.cut_names.append('%s_%s'%(self.prefix,cs))
        self.reset()
        
        file_out.cd()
        self.hEvents = ROOT.TH1I('hEvents_%s'%(self.name), '', 1+len(self.cut_names), -0.5, 0.5+len(self.cut_names))
        self.hEvents.GetXaxis().SetTitle('')
        self.hEvents.GetYaxis().SetTitle('nEvents')
        self.hEvents.GetYaxis().SetTitleOffset(1.25)
        for icn in range(0,len(self.cut_suffixes)):
            self.hEvents.GetXaxis().SetBinLabel(icn+2,self.cut_suffixes[icn])
        self.hEvents.GetXaxis().SetBinLabel(1,'raw')
        self.hEvents.SetFillColor(ROOT.kRed)
        self.hEvents.SetLineColor(ROOT.kRed)
        self.hEvents.SetMarkerColor(ROOT.kRed)
        self.hEvents.SetMarkerStyle(20)
        self.hEvents.SetMinimum(0)
        
        self.nFail = 0
    
    def analyse_events(self, events, sname):
        # Make histograms
        h_raw = {}
        h_cum = {}
        h_cut = {}
        h_NM1 = {}
        file_out.cd()
        for cs in cut_suffixes[self.suffix]:
            h_raw[cs] = {}
            h_cum[cs] = {}
            h_cut[cs] = {}
            h_NM1[cs] = {}
            for rname in regions:
                h_raw[cs][rname] = vars[cs].histogram.Clone('hraw_%s_%s_%s_%s'%(self.name,cs,sname,rname))
                h_cum[cs][rname] = vars[cs].histogram.Clone('hcum_%s_%s_%s_%s'%(self.name,cs,sname,rname))
                h_cut[cs][rname] = vars[cs].histogram.Clone('hcut_%s_%s_%s_%s'%(self.name,cs,sname,rname))
                h_NM1[cs][rname] = vars[cs].histogram.Clone('hNM1_%s_%s_%s_%s'%(self.name,cs,sname,rname))
                h_cum[cs][rname].SetFillColor(ROOT.kYellow)
        
        self.reset()
        for ev in events:
            ev_success            = {}
            ev_success_cumulative = {}
            ev_success_NM1        = {}
            self.nEvents['raw']    += 1
            self.nElectrons['raw'] += len(ev)
            for cn in self.cut_names:
                ev_success[cn]            = 0
                ev_success_NM1[cn]        = 0
                ev_success_cumulative[cn] = 0
            for i_cn in range(0,len(self.cut_names)):
                cn = self.cut_names[i_cn]
                hKey = cut_suffixes[self.suffix][i_cn]
                for el in ev:
                    if el.isValid == False:
                        continue
                    if not el.cutflow_bitmaps:
                        continue
                    
                    bitmap = el.cutflow_bitmaps[self.name]
                    b_thisCut           = int(math.pow(2,i_cn+1)  )
                    b_thisCutCumulative = int(math.pow(2,i_cn+1)-1)
                    b_thisCutNM1        = int(math.pow(2,len(self.cut_names)+1)-math.pow(2,i_cn+1)-1)
                    
                    value = el.cutflows[self.name][i_cn+1].value
                    r = el.region
                    h_raw[hKey][r].Fill(value)
                    if (bitmap & b_thisCut) == b_thisCut:
                        ev_success[cn]            += 1
                        self.nElectrons[cn] += 1
                        h_cut[hKey][r].Fill(value)
                    if (bitmap & b_thisCutCumulative) == b_thisCutCumulative:
                        ev_success_cumulative[cn] += 1
                        h_cum[hKey][r].Fill(value)
                        self.cumulative_nElectrons[cn] += 1
                    if (bitmap & b_thisCutNM1) == b_thisCutNM1:
                        ev_success_NM1     [cn] += 1
                        self.NM1_nElectrons[cn] += 1
                        h_NM1[hKey][r].Fill(value)
            for cn in self.cut_names:
                if ev_success[cn] >= 2:
                    self.nEvents[cn] += 1
                if ev_success_cumulative[cn] >= 2:
                    self.cumulative_nEvents[cn] += 1
                if ev_success_NM1[cn] >= 2:
                    self.NM1_nEvents[cn] += 1
        
        # Save histograms to file
        file_out.cd()
        for hName in h_raw:
            for rname in regions:
                h_raw[hName][rname].Write()
                h_cut[hName][rname].Write()
                h_cum[hName][rname].Write()
                h_NM1[hName][rname].Write()
   
    def print_results(self):
        print '%40s  %15s  %15s  %15s  %15s  %15s  %15s'%('', 'nEvents', 'cum nEvents', 'N-1 events', 'nElectrons', 'cum nElectrons', 'N-1 electrons')
        print '%40s  %15s  %15d  %15s  %15s  %15d  %15s'%('Raw', '', self.nEvents['raw'], '', '', self.nElectrons['raw'], '')
        for cn in self.cut_names:
            print '%40s  %15d  %15d  %15d  %15d  %15d  %15d'%(cn, self.nEvents[cn], self.cumulative_nEvents[cn], self.NM1_nEvents[cn], self.nElectrons[cn], self.cumulative_nElectrons[cn], self.NM1_nElectrons[cn])
    def make_histogram(self, sname):
        h = self.hEvents.Clone('hEvents_%s_%s'%(self.name,sname))
        for icn in range(0,len(self.cut_names)):
            h.SetBinContent(icn+2,self.cumulative_nEvents[self.cut_names[icn]])
        h.SetBinContent(1,self.nEvents['raw'])
        return h
    def reset(self):
        for cn in self.cut_names:
            self.nEvents[cn]               = 0
            self.NM1_nEvents[cn]           = 0
            self.cumulative_nEvents[cn]    = 0
            self.nElectrons[cn]            = 0
            self.NM1_nElectrons[cn]        = 0
            self.cumulative_nElectrons[cn] = 0
        self.nEvents['raw']    = 0
        self.nElectrons['raw'] = 0
        
cutflows = []
for cp in cutflow_prefixes:
    for cs in collection_suffixes:
        if 'total' not in cs:
            continue
        cutflows.append(cutflow(cp, cs))

sname = 'DYEE'
sname = 'ZprimeToEE_M5000_v5'
file = ROOT.TFile('ntuples/outfile_%s.root'%sname,'READ')
#file = ROOT.TFile('ntuples/20141002/outfile_%s_skimmed_slimmed.root'%sname,'READ')
tree = file.Get('IIHEAnalysis')
tree.SetBranchStatus('*'     ,0)
tree.SetBranchStatus('mc_*'  ,1)
tree.SetBranchStatus('gsf_*' ,1)
tree.SetBranchStatus('HEEP_*',1)

events = []
nEntries = tree.GetEntries()
nEntries = 10000
for i in range(0, nEntries):
    if i%1000==0:
        print '%8d / %8d' %(i,nEntries)
    tree.GetEntry(i)
        
    electrons = []
    for i_el in range(0, tree.gsf_n):
        el = electron_object(cutflows, tree, i_el)
        if el.isValid:
            electrons.append(el)
    events.append(electrons)
    
for c in cutflows:
    cname = c.name
    if 'total' not in cname:
        continue
    c.analyse_events(events, sname)
    file_out.cd()
    h = c.make_histogram(sname)
    c.print_results()
    h.Draw('e3')
    canvas.Print('plots/%s_%s.eps'%(cname, sname))
    canvas.Print('plots/%s_%s.png'%(cname, sname))
    
    file_out.cd()
    h.Write()
    print

file_out.Write()
