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

file_in = ROOT.TFile('histograms.root','READ')

sname = 'ZprimeToEE_M5000_v5'
#sname = 'DYEE'
sample_names = []
sample_names.append('ZprimeToEE_M5000_v5')
#sample_names.append('DYEE')

cutflow_names = []
cutflow_names.append('HEEP_cutflow41_total'     )
cutflow_names.append('HEEP_cutflow50_50ns_total')
cutflow_names.append('HEEP_cutflow50_25ns_total')

colors = {}
colors['HEEP_cutflow41_total'     ] = ROOT.kRed
colors['HEEP_cutflow50_50ns_total'] = ROOT.kBlue
colors['HEEP_cutflow50_25ns_total'] = ROOT.kGreen


h_cutflows = {}
for cname in cutflow_names:
    h_cutflows[cname] = {}
    for sname in sample_names:
        h_cutflows[cname][sname] = file_in.Get('hEvents_%s_%s'%(cname,sname))

first = True
legx = 0.2
legy = 0.2
legend = ROOT.TLegend(legx, legy, legx+0.4, legy+0.2)
legend.SetFillColor(0)
legend.SetShadowColor(0)
legend.SetBorderSize(0)
for cname in cutflow_names:
    draw_options = 'e3' if first else 'e3:sames'
    h = h_cutflows[cname][sname]
    denom = 1.0*h.GetBinContent(1)
    for bin in range(1, h.GetNbinsX()+1):
        content = h.GetBinContent(bin)
        eff = 1.0*content/denom
        err = math.sqrt(eff*(1-eff)/denom)
        #print bin , content , denom , eff , err
        h.SetBinContent(bin, eff*100)
        h.SetBinError  (bin, err*100)
    h.GetYaxis().SetTitle('cumulative effiency (%)')
    h.SetFillColor(colors[cname])
    h.SetLineColor(colors[cname])
    h.SetMarkerColor(colors[cname])
    #h.SetMaximum(1000)
    legend.AddEntry(h, cname, 'f')
    
    h.Draw(draw_options)
    first = False
legend.Draw()
canvas.Print('plots/hCutflow_%s.eps'%sname)

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

canvas_2plots = ROOT.TCanvas('canvas_2plots', '', 100, 100, 1200, 600)
canvas_2plots.Divide(2,1)
canvas_2plots.GetPad(1).SetGridx()
canvas_2plots.GetPad(1).SetGridy()
canvas_2plots.GetPad(2).SetGridx()
canvas_2plots.GetPad(2).SetGridy()



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

regions = ['barrel','endcap']

for cname in cutflow_names:
    for vname in var_names:
        for rname in regions:
            hRaw = file_in.Get('hraw_%s_%s_%s_%s'%(cname, vname, sname, rname))
            hCum = file_in.Get('hcum_%s_%s_%s_%s'%(cname, vname, sname, rname))
            hNM1 = file_in.Get('hNM1_%s_%s_%s_%s'%(cname, vname, sname, rname))
            
            legend = ROOT.TLegend(0.12,0.85,0.88,0.78)
            legend.SetFillColor(0)
            legend.SetBorderSize(0)
            legend.SetShadowColor(0)
            legend.SetNColumns(2)
            legend.AddEntry(hRaw, 'Raw distribution'       , 'f')
            legend.AddEntry(hCum, 'Cumulative distribution', 'f')
        
            canvas_2plots.GetPad(1).SetLogy(var_log[vname])
            canvas_2plots.GetPad(2).SetLogy(var_log[vname])
        
            hRaw.GetYaxis().SetTitleOffset(1.5)
            hCum.GetYaxis().SetTitleOffset(1.5)
            hNM1.GetYaxis().SetTitleOffset(1.5)
        
            multiplier = 10 if var_log[vname] else 1.25
            hRaw.SetMaximum(multiplier*hRaw.GetMaximum())
            hCum.SetMaximum(multiplier*hCum.GetMaximum())
            hNM1.SetMaximum(multiplier*hNM1.GetMaximum())
        
            canvas_2plots.cd(1)
            hRaw.Draw()
            hCum.Draw('sames')
            hRaw.Draw('sames:axis')
            legend.Draw()
            CMS_label.Draw()
        
            canvas_2plots.cd(2)
            hNM1.Draw()
            NM1_label.Draw()
            CMS_label.Draw()
        
            canvas_2plots.Print('plots/vars/%s_%s_%s_%s.eps'%(cname, vname, sname, rname))

