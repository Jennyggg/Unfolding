import CMS_lumi
import ROOT as rt; rt.PyConfig.IgnoreCommandLineOptions = True
import array
import numpy as np

histLineColor = rt.kBlack
markerSize  = 1.0
histFillColor =['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#a65628','#f781bf','#999999','#ffff33']


#CMS Style
iPeriod = 0
iPos = 11
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Preliminary"
CMS_lumi.lumi_sqrtS = "35.9 fb^{-1} (13 TeV)" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
if( iPos==0 ): CMS_lumi.relPosX = 0.12
H_ref = 600;
W_ref = 540;
x1_l = 0.95
y1_l = 0.90
dx_l = 0.30
dy_l = 0.28
x0_l = x1_l-dx_l
y0_l = y1_l-dy_l
ar_l = dy_l/dx_l
W = W_ref
H  = H_ref
# references for T, B, L, R
T = 0.08*H_ref
B = 0.12*H_ref
L = 0.14*W_ref
R = 0.04*W_ref

H_long = 600;
W_long = 700;

gap_ = 1.0/6
bwx_ = 0.14/2
bwy_ = gap_/1.5

x_l = [1.2*bwx_]

y_l = [1-gap_]
ex_l = [0]
ey_l = [0.04/ar_l]


x_l = array.array("f",x_l)
ex_l = array.array("f",ex_l)
y_l = array.array("f",y_l)
ey_l = array.array("f",ey_l)
xx_ = x_l[0]



#########################################################
#Plot specific quantities
#Ratio specific
RTSX = 0.156 #ratio title size
RTOX = 1.1 #ratio title offset
RLSX = 0.161 #ratio label size
RTSY = 0.166 #ratio title size
RTOY = 0.3 #ratio title offset
RLSY = 0.138 #ratio label size
RMIN = 0.7
RMAX = 1.3


#Second pad quantities
P2RM = 0.04 #pad right margin
P2LM = 0.12 #pad left margin
P2TM = 0.05  #pad top margin
P2BM = 0.42  #pad bottom margin
#Legend specific
xgap = -0.46
n_ = 6 #Number of items on legend
gap_ = 1./(n_+1)
#gap_ = 0.235
#ygap = -0.45
xshiftm = -0.2
xshiftp = -0.0
yshiftm =+0.0
yshiftp = 0.00
legendsize = 0.09
#Frame specific
FLS = 0.04
FTS = 0.055
FTO = 0.98
FRML = 2.0 #linear maximum frame range multiplyier
FRMLOG = 5e4 #log maximum frame range multiplyier
ranges = 1000
bins = 50


def SetupBox(box_,yy_,fill = rt.kBlack):

    box_.SetLineStyle( rt.kSolid )
    box_.SetLineWidth( 1 )
    box_.SetLineColor( rt.kBlack )
    box_.SetFillColor(fill)


def SetupCanvas(c,logy):
    c.SetFillColor(0)
    c.SetBorderMode(0)
    c.SetFrameFillStyle(0)
    c.SetFrameBorderMode(0)
    c.SetLeftMargin( L/W_ref)
    c.SetRightMargin( R/W_ref)
    c.SetTopMargin( T/H )
    c.SetBottomMargin( B/H )
    c.SetTickx(0)
    c.SetTicky(0)
    if logy: c.SetLogy()


def AddOverflow(h):
    b0 = h.GetBinContent(0)
    e0 = h.GetBinError(0)
    nb = h.GetNbinsX()
    bn = h.GetBinContent(nb + 1)
    en = h.GetBinError(nb + 1)

    h.SetBinContent(0, 0)
    h.SetBinContent(nb+1, 0)
    h.SetBinError(0, 0)
    h.SetBinError(nb+1, 0)

    h.SetBinContent(1, h.GetBinContent(1) + b0)
    h.SetBinError(1, (h.GetBinError(1)**2 + e0**2)**0.5 )

    h.SetBinContent(nb, h.GetBinContent(nb) + bn)
    h.SetBinError(nb, (h.GetBinError(nb)**2 + en**2)**0.5 )

def AddOverflow2D(h):
    nbx = h.GetNbinsX()
    nby = h.GetNbinsY()
    for ix in range(nbx):
        b0x=h.GetBinContent(ix,0)
        e0x=h.GetBinError(ix,0)
        h.SetBinContent(ix,0, 0)
        h.SetBinError(ix,0, 0)
        h.SetBinContent(ix,1, h.GetBinContent(ix,1) + b0x)
        h.SetBinError(ix,1, (h.GetBinError(ix,1)**2 + e0x**2)**0.5 )


        bx = h.GetBinContent(ix,nby+1)
        ex = h.GetBinContent(ix,nby+1)
        h.SetBinContent(ix,nby+1, 0)
        h.SetBinError(ix,nby+1, 0)
        h.SetBinContent(ix,nby, h.GetBinContent(ix,nby) + bx)
        h.SetBinError(ix,nby, (h.GetBinError(ix,nby)**2 + ex**2)**0.5 )

    for iy in range(nby):
        b0y=h.GetBinContent(0,iy+1)
        e0y=h.GetBinError(0,iy+1)
        h.SetBinContent(0,iy+1,0)
        h.SetBinError(0,iy+1,0)
        h.SetBinContent(1,iy+1, h.GetBinContent(1,iy+1) + b0y)
        h.SetBinError(1,iy+1, (h.GetBinError(1,iy+1)**2 + e0y**2)**0.5 )


        by = h.GetBinContent(nbx+1,iy+1)
        ey = h.GetBinContent(nbx+1,iy+1)
        h.SetBinContent(nbx+1,ix+1, 0)
        h.SetBinError(nbx+1,iy+1, 0)
        h.SetBinContent(nbx,iy+1, h.GetBinContent(nbx,iy+1) + by)
        h.SetBinError(nbx,iy+1, (h.GetBinError(nbx,iy+1)**2 + ey**2)**0.5 )

def create_paves(lumi, label, CMSposX=0.11, CMSposY=0.9, prelimPosX=0.11, prelimPosY=0.85, lumiPosX=0.95, lumiPosY=0.951, alignRight=False, CMSsize=0.75*0.08, prelimSize=0.75*0.08*0.76, lumiSize=0.6*0.08):

    #pt_lumi = ROOT.TPaveText(xhi-0.25, ylo, xhi, yhi,"brNDC")
    pt_lumi = rt.TPaveText(lumiPosX-0.25, lumiPosY, lumiPosX, 1.0,"brNDC")
    pt_lumi.SetFillStyle(0)
    pt_lumi.SetBorderSize(0)
    pt_lumi.SetFillColor(0)
    pt_lumi.SetTextFont(42)
    pt_lumi.SetTextSize(lumiSize)
    pt_lumi.SetTextAlign(31) #left=10, bottom=1, centre=2
 #input lumi in fb^-1
    if np.log10(lumi)>=0 and np.log10(lumi)<3:
      lumi_print = lumi
      lumi_unit = 'fb^{-1}'
    elif np.log10(lumi)<6 and np.log10(lumi)>=3:
      lumi_print = lumi*0.001
      lumi_unit = 'ab^{-1}'
    elif np.log10(lumi)<0 and np.log10(lumi)>=-3:
      lumi_print = lumi*1000
      lumi_unit = 'pb^{-1}'
    elif np.log10(lumi)<-3 and np.log10(lumi)>=-6:
      lumi_print = lumi*1e6
      lumi_unit = 'nb^{-1}'
    elif np.log10(lumi)<=-6 and np.log10(lumi)>=-9:
      lumi_print = lumi*1e9
      lumi_unit = '#mub^{-1}'
    elif np.log10(lumi)<=-9 and np.log10(lumi)>=-12:
      lumi_print = lumi*1e12
      lumi_unit = 'mb^{-1}'
    elif np.log10(lumi)<=-12 and np.log10(lumi)>=-15:
      lumi_print = lumi*1e15
      lumi_unit = 'b^{-1}'


    pt_lumi.AddText( "{0:.1f}".format(lumi_print)+" "+lumi_unit+" (13 TeV)" )

    # if CMSpos == 0: #outside frame
    #     pt_CMS = rt.TPaveText(xlo, ylo, xlo+0.1, yhi,"brNDC")
    # elif CMSpos == 1: #left
    #     pt_CMS = rt.TPaveText(xlo+0.04, ylo-0.09, xlo+0.14, ylo-0.04,"brNDC")
    # elif CMSpos == 2: #center
    #     pt_CMS = rt.TPaveText(xlo+0.4, ylo-0.09, xlo+0.5, ylo-0.04,"brNDC")
    # elif CMSpos == 3: #right
    #     pt_CMS = rt.TPaveText(xhi-0.2, ylo-0.09, xhi-0.1, ylo-0.04,"brNDC")
    if alignRight:
        pt_CMS = rt.TPaveText(CMSposX-0.1, CMSposY, CMSposX, CMSposY+0.05,"brNDC")
    else:
        pt_CMS = rt.TPaveText(CMSposX, CMSposY, CMSposX+0.1, CMSposY+0.05,"brNDC")
    pt_CMS.SetFillStyle(0)
    pt_CMS.SetBorderSize(0)
    pt_CMS.SetFillColor(0)
    pt_CMS.SetTextFont(61)
    pt_CMS.SetTextSize(CMSsize)
    #pt_CMS.SetTextAlign(31 if CMSpos==3 else 11)
    pt_CMS.SetTextAlign(31 if alignRight else 11 )
    pt_CMS.AddText("CMS")

    # if PrelimPos == 0: #outside frame
    #     pt_prelim = rt.TPaveText(xlo+0.09, ylo, xlo+0.3, yhi,"brNDC")
    # elif PrelimPos == 1: #left beside CMS
    #     pt_prelim = rt.TPaveText(xlo+0.13, ylo-0.09, xlo+0.34, ylo-0.04,"brNDC")
    # elif PrelimPos == 2: #left under CMS
    #     pt_prelim = rt.TPaveText(xlo+0.04, ylo-0.15, xlo+0.14, ylo-0.10,"brNDC")
    # elif PrelimPos == 3: #right under CMS
    #     pt_prelim = rt.TPaveText(xhi-0.2, ylo-0.15, xhi-0.1, ylo-0.10,"brNDC")
    if alignRight:
        pt_prelim = rt.TPaveText(prelimPosX-0.2, prelimPosY, prelimPosX, prelimPosY+0.05,"brNDC")
    else:
        pt_prelim = rt.TPaveText(prelimPosX, prelimPosY, prelimPosX+0.2, prelimPosY+0.05,"brNDC")
    pt_prelim.SetFillStyle(0)
    pt_prelim.SetBorderSize(0)
    pt_prelim.SetFillColor(0)
    pt_prelim.SetTextFont(52)
    pt_prelim.SetTextSize(prelimSize)
    #pt_prelim.SetTextAlign(31 if PrelimPos==3 else 11)
    pt_prelim.SetTextAlign(31 if alignRight else 11 )
    if label == "SimPAS":
        pt_prelim.AddText("Simulation Preliminary")
    elif label == "DataPAS":
        pt_prelim.AddText("Preliminary")
    elif label == "WP":
        pt_prelim.AddText("Private work")
    elif label == "Sim":
        pt_prelim.AddText("Simulation")
    elif label == "Data":
        pt_prelim.AddText("")
    elif label == "SimSupp":
        pt_prelim.AddText("Simulation Supplementary")
    elif label == "DataSupp":
        pt_prelim.AddText("Supplementary")
    elif label == "Int":
        pt_prelim.AddText("Internal")

    return {"lumi":pt_lumi, "CMS":pt_CMS, "label":pt_prelim}

