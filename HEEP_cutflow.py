##########################################################################################
#                                  HEEP cutflow plotter                                  #
##########################################################################################
# (c) 2015 Aidan Randle-Conde (ULB)                                                      #
# Contact: aidan.randleconde@gmail.com                                                   #
# Github:  https://github.com/ULBHEEPAnalyses/HEEPCutFlow                                #
##########################################################################################

import math

##########################################################################################
#                                  Settings for the job                                  #
##########################################################################################
sname = 'ZprimeToEE_M5000_v5'

##########################################################################################
#                             Import ROOT and apply settings                             #
##########################################################################################
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

##########################################################################################
#                           Create canvas, create output file                            #
##########################################################################################
canvas = ROOT.TCanvas('canvas', '', 100, 100, 1200, 800)
canvas.SetGridx()
canvas.SetGridy()

file_out = ROOT.TFile('histograms.root','RECREATE')

##########################################################################################
#                                Create variable objects                                 #
##########################################################################################
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

def add_variable(vname, nBins, lower, upper, xaxis, unit):
    var_names.append(vname)
    vars[vname] = var_object(vname, nBins, lower, upper, xaxis, unit)

add_variable('Et'             ,  100,     0, 3000, 'E_{T}(e)'                , 'GeV')
add_variable('eta'            ,  100,    -3,    3, '#eta_{T}(e)'             , ''   )
add_variable('EcalDriven'     ,    3,  -1.5,  1.5, 'isEcalDriven(e)'         , ''   )
add_variable('dEtaIn'         ,  100, -0.05, 0.05, '#Delta#eta_{in}(e)'      , ''   )
add_variable('dPhiIn'         ,  100,  -0.1,  0.1, '#Delta#phi_{in}(e)'      , ''   )
add_variable('HOverE'         ,  100,     0,    1, 'H/E(e)'                  , ''   )
add_variable('SigmaIetaIeta'  ,  100,     0, 0.05, '#sigma_{i#etai#eta}(e)'  , ''   )
add_variable('E1x5OverE5x5'   ,  100,     0,    1, 'E_{1x5}/E_{5x5}(e)'      , ''   )
add_variable('E2x5OverE5x5'   ,  100,   0.5,    1, 'E_{2x5}/E_{5x5}(e)'      , ''   )
add_variable('missingHits'    ,   11,  -0.5, 10.5, 'missing hits (e)'        , ''   )
add_variable('dxyFirstPV'     ,  100,  -0.1,  0.1, 'd_{xy}^{PV}(e)'          , 'mm' )
add_variable('isolEMHadDepth1',  100,     0,   50, 'isolation_{HadDepth1}(e)', 'GeV')
add_variable('IsolPtTrks'     ,  100,     0,   10, 'isolation_{Tracks}(e)'   , ''   )

##########################################################################################
#                                 Assemble the cutflows                                  #
# This is somewhat arbitrary.  The HEEP group does not define an order of cuts.          #
##########################################################################################
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

##########################################################################################
#                              Class for storing electrons                               #
##########################################################################################
regions = ['barrel','endcap']

class cutflow_instance_object:
    def __init__(self, name, value, success):
        self.name    = name
        self.value   = value
        self.success = success

class electron_object:
    def __init__(self, cutflows_in, tree, i_el):
        # Quick sanity check
        self.isValid = i_el < len(tree.HEEP_cutflow41_total)
        if self.isValid == False:
            return
        
        # Everthing is okay, so store p4 and information about detector region
        self.pt     = tree.gsf_pt[i_el]
        self.eta    = tree.gsf_eta[i_el]
        self.phi    = tree.gsf_phi[i_el]
        self.energy = tree.gsf_energy[i_el]
        self.region = 'barrel' if abs(self.eta) < 1.4442 else 'endcap'
        self.p4 = ROOT.TLorentzVector()
        self.p4.SetPtEtaPhiE(self.pt, self.eta, self.phi, self.energy)
        
        # Match against the truth record to ensure E(reco)/E(gen) > 0.9
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
        
        # For each electron, analyse the various cutflows
        # Results are stored as bitmaps which are parsed later on
        # So a given cutflow might look like "110100110"
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

##########################################################################################
#                              Class for storing cutflows                                #
# This class does all the hard work.  For each cut it stores:                            #
#   - the number of events that pass that cut                                            #
#   - the number of electrons that pass that cut                                         #
# For both events and electrons the numbers are stored in terms of:                      #
#   - passing that particular cut (cut n)                                                #
#   - passing that particular cut and all before it (cumulative n)                       #
#   - passing all but that particular cut (N-1 n)                                        #
# The class also has a histogram for showing the distributions of each variables for     #
# four different scenarios:                                                              #
#   - no selection applied (raw distribution)                                            #
#   - passing that particular cut (cut distribution)                                     #
#   - passing that particular cut and all before it (cumulative distribution)            #
#   - passing all but that particular cut (N-1 distribution)                             #
##########################################################################################
class cutflow:
    def __init__(self, prefix, suffix):
        self.prefix = prefix
        self.suffix = suffix
        self.name   = '%s_%s'%(self.prefix,self.suffix)
        self.cut_names = []
        
        # Store information about how many events/electrons pass each cut
        self.nEvents               = {}
        self.NM1_nEvents           = {}
        self.cumulative_nEvents    = {}
        self.nElectrons            = {}
        self.NM1_nElectrons        = {}
        self.cumulative_nElectrons = {}
        
        # Get the cut names.  Somewhat clunky way to do this.
        self.cut_suffixes = cut_suffixes[self.suffix]
        for cs in cut_suffixes[self.suffix]:
            self.cut_names.append('%s_%s'%(self.prefix,cs))
        self.reset()
        self.nFail = 0
        
        # Create histogram to store events, save this to file
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
    
    def analyse_events(self, events, sname):
        # Make histograms, each with a unique name
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
        
        # Reset values and loop over events
        self.reset()
        for ev in events:
            ev_success            = {}
            ev_success_cumulative = {}
            ev_success_NM1        = {}
            
            # Increment the raw values. (Recall that an "event" is a list of "electron" classes.)
            self.   nEvents['raw'] += 1
            self.nElectrons['raw'] += len(ev)
            
            # Reset the event success flag for this event for each cut
            for cn in self.cut_names:
                ev_success[cn]            = 0
                ev_success_NM1[cn]        = 0
                ev_success_cumulative[cn] = 0
            
            # Now loop over the cuts to get their flag for this event
            for i_cn in range(0,len(self.cut_names)):
                cn = self.cut_names[i_cn]              # Cut name
                hKey = cut_suffixes[self.suffix][i_cn] # array key for histogram
                
                # Loop over the electrons in this event
                for el in ev:
                    # Sanity checks to make sure the electron is "good" and has cutflow bitmaps stored
                    if el.isValid == False:
                        continue
                    if not el.cutflow_bitmaps:
                        continue
                    
                    # Now parse the bitmaps.  This is a bit messy.
                    bitmap = el.cutflow_bitmaps[self.name]
                    b_thisCut           = int(math.pow(2,i_cn+1)  )
                    b_thisCutCumulative = int(math.pow(2,i_cn+1)-1)
                    b_thisCutNM1        = int(math.pow(2,len(self.cut_names)+1)-math.pow(2,i_cn+1)-1)
                    
                    # Get the value of the variable for this electron
                    value = el.cutflows[self.name][i_cn+1].value
                    r = el.region
                    
                    # Now fill the histograms
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
            
            # Fill in the cutflow information for the event.
            # An event passes if at least two electrons pass.
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
        # Print results to screen in a nice table
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
            self.              nEvents[cn] = 0
            self.          NM1_nEvents[cn] = 0
            self.   cumulative_nEvents[cn] = 0
            self.           nElectrons[cn] = 0
            self.       NM1_nElectrons[cn] = 0
            self.cumulative_nElectrons[cn] = 0
        self.   nEvents['raw'] = 0
        self.nElectrons['raw'] = 0

##########################################################################################
#                                 Assemble the cutflows                                  #
# This is somewhat arbitrary.  The HEEP group does not define an order of cuts.          #
##########################################################################################
cutflow_prefixes = []
cutflow_prefixes.append('HEEP_cutflow41'     )
cutflow_prefixes.append('HEEP_cutflow50_25ns')
cutflow_prefixes.append('HEEP_cutflow50_50ns')

cutflows = []
for cp in cutflow_prefixes:
    for cs in collection_suffixes:
        if 'total' not in cs:
            continue
        cutflows.append(cutflow(cp, cs))

##########################################################################################
#                    Read the input file, turn off unwanted branches                     #
##########################################################################################
file = ROOT.TFile('ntuples/outfile_%s.root'%sname,'READ')
tree = file.Get('IIHEAnalysis')
tree.SetBranchStatus('*'     ,0)
tree.SetBranchStatus('mc_*'  ,1)
tree.SetBranchStatus('gsf_*' ,1)
tree.SetBranchStatus('HEEP_*',1)

##########################################################################################
#                      Run over the events to create the electrons                       #
##########################################################################################
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

##########################################################################################
#      Populated the cutflows with the events, save to file, and print raw results       #
##########################################################################################
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
