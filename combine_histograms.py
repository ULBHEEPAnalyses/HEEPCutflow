##########################################################################################
#                                  HEEP cutflow plotter                                  #
##########################################################################################
# (c) 2015 Aidan Randle-Conde (ULB)                                                      #
# Contact: aidan.randleconde@gmail.com                                                   #
# Github:  https://github.com/ULBHEEPAnalyses/HEEPCutFlow                                #
##########################################################################################

import math
import ROOT

##########################################################################################
#                                  Settings for the job                                  #
##########################################################################################
sname = 'ZprimeToEE_M5000_v5'

colors = {}
colors['HEEP_cutflow41_total'     ] = ROOT.kRed
colors['HEEP_cutflow50_50ns_total'] = ROOT.kBlue
colors['HEEP_cutflow50_25ns_total'] = ROOT.kGreen

##########################################################################################
#                             Import ROOT and apply settings                             #
##########################################################################################
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
#                                 Create canvas, labels                                  #
##########################################################################################
canvas = ROOT.TCanvas('canvas', '', 100, 100, 1200, 800)
canvas.SetGridx()
canvas.SetGridy()

canvas_plots = ROOT.TCanvas('canvas_plots', '', 100, 100, 1200, 600)
canvas_plots.Divide(2,1)
canvas_plots.GetPad(1).SetGridx()
canvas_plots.GetPad(1).SetGridy()
canvas_plots.GetPad(2).SetGridx()
canvas_plots.GetPad(2).SetGridy()

NM1_label = ROOT.TLatex(0.2,0.8,'(N-1) distribution')
NM1_label.SetNDC()

CMS_label_texts = {}
CMS_label_texts['normal'        ] = 'CMS'
CMS_label_texts['internal'      ] = 'CMS internal'
CMS_label_texts['workInProgress'] = 'CMS work in progress'
CMS_labels = {}
for t in CMS_label_texts:
    CMS_labels[t] = ROOT.TLatex(0.65, 0.945, CMS_label_texts[t])
    CMS_labels[t].SetNDC()
CMS_label = CMS_labels['internal']

##########################################################################################
#                             Get input file, set up cutflows                            #
##########################################################################################
file_in = ROOT.TFile('histograms.root','READ')

sample_names = []
sample_names.append(sname)

cutflow_names = []
cutflow_names.append('HEEP_cutflow41_total'     )
cutflow_names.append('HEEP_cutflow50_50ns_total')
cutflow_names.append('HEEP_cutflow50_25ns_total')

h_cutflows = {}
for cname in cutflow_names:
    h_cutflows[cname] = {}
    for sname in sample_names:
        h_cutflows[cname][sname] = file_in.Get('hEvents_%s_%s'%(cname,sname))

first = True

##########################################################################################
#                                     Create legend                                      #
##########################################################################################
legx = 0.2
legy = 0.2
legend_cutflow = ROOT.TLegend(legx, legy, legx+0.4, legy+0.2)
legend_cutflow.SetFillColor(0)
legend_cutflow.SetShadowColor(0)
legend_cutflow.SetBorderSize(0)

##########################################################################################
#                                   Draw everything                                      #
##########################################################################################
canvas.cd()
for cname in cutflow_names:
    draw_options = 'e3' if first else 'e3:sames'
    h = h_cutflows[cname][sname]
    denom = 1.0*h.GetBinContent(1)
    for bin in range(1, h.GetNbinsX()+1):
        content = h.GetBinContent(bin)
        eff = 1.0*content/denom
        err = math.sqrt(eff*(1-eff)/denom)
        h.SetBinContent(bin, eff*100)
        h.SetBinError  (bin, err*100)
    h.GetYaxis().SetTitle('cumulative effiency (%)')
    h.SetFillColor(colors[cname])
    h.SetLineColor(colors[cname])
    h.SetMarkerColor(colors[cname])
    legend_cutflow.AddEntry(h, cname, 'f')
    
    h.Draw(draw_options)
    first = False
legend_cutflow.Draw()
canvas.Print('plots/hCutflow_%s.eps'%sname)

##########################################################################################
#                              Plot variable distributions                               #
##########################################################################################
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

var_log = {}
var_log_names = ['Et','dEtaIn','dPhiIn','HOverE','SigmaIetaIeta','dxyFirstPV','IsolPtTrks']
for vname in var_names:
    var_log[vname] = False
for vname in var_log_names:
    var_log[vname] = True

regions = ['barrel','endcap']

for cname in cutflow_names:
    for vname in var_names:
        for rname in regions:
            hRaw = file_in.Get('hraw_%s_%s_%s_%s'%(cname, vname, sname, rname))
            hCum = file_in.Get('hcum_%s_%s_%s_%s'%(cname, vname, sname, rname))
            hNM1 = file_in.Get('hNM1_%s_%s_%s_%s'%(cname, vname, sname, rname))
            
            legend_var = ROOT.TLegend(0.12,0.85,0.88,0.78)
            legend_var.SetFillColor(0)
            legend_var.SetBorderSize(0)
            legend_var.SetShadowColor(0)
            legend_var.SetNColumns(2)
            legend_var.AddEntry(hRaw, 'Raw distribution'       , 'f')
            legend_var.AddEntry(hCum, 'Cumulative distribution', 'f')
            
            canvas_plots.GetPad(1).SetLogy(var_log[vname])
            canvas_plots.GetPad(2).SetLogy(var_log[vname])
        
            hRaw.GetYaxis().SetTitleOffset(1.5)
            hCum.GetYaxis().SetTitleOffset(1.5)
            hNM1.GetYaxis().SetTitleOffset(1.5)
        
            multiplier = 10 if var_log[vname] else 1.25
            hRaw.SetMaximum(multiplier*hRaw.GetMaximum())
            hCum.SetMaximum(multiplier*hCum.GetMaximum())
            hNM1.SetMaximum(multiplier*hNM1.GetMaximum())
            
            # Draw everything
            # Raw and cumulative first
            canvas_plots.cd(1)
            hRaw.Draw()
            hCum.Draw('sames')
            hRaw.Draw('sames:axis')
            legend_var.Draw()
            CMS_label.Draw()
            
            # Then N-1
            canvas_plots.cd(2)
            hNM1.Draw()
            NM1_label.Draw()
            CMS_label.Draw()
        
            canvas_plots.Print('plots/vars/%s_%s_%s_%s.eps'%(cname, vname, sname, rname))

