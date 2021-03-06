#mac user that experience some problems can try export VERSIONER_PYTHON_PREFER_32_BIT=yes BEFORE THE SCRIPT!!!!

import wx
import os
import sys
import subprocess
import ctypes
# create simple class: Just a window. If you come from VB, a frame here is like VB window
#input=''
#probe=''

input_sad=''

class MainFrame(wx.Frame):
    
    
    
    def __init__(self):
        wx.Frame.__init__(self, None,wx.ID_ANY, title='PIPE-ChipSAD', size=(900, 630))
        # self (which is now Text editor ) is now a Window with all attribute of a frame
        # Your widgets like buttons, text controls etc will be added here before showing the frame
        # as example I will add a panel
        
        
        self.background = wx.Panel(self)
        #self.SetBackgroundColour((134, 236, 789))
        #self.inputArea = wx.TextCtrl(self.background)
        
        img=wx.Image( "./imgs_chipsad_gui/chipSAD_logo.jpg",  wx.BITMAP_TYPE_ANY).ConvertToBitmap()
        self.bitmap1 = wx.StaticBitmap(self.background, -1, img, (0, 0))#, size=(400, 150))
        
        ##img1=wx.Image( "./imgs_chipsad_gui/punto_esclam_rosso.jpeg",  wx.BITMAP_TYPE_ANY).ConvertToBitmap()
        ##self.bitmap2 = wx.StaticBitmap(self.background, -1, img1, (500, 50))#, size=(400, 150))
        
        wx.StaticText(self.background, -1, 'Microarray design options:', (500, 60) )
        
        wx.StaticText(self.background, 1, 'Choose the signal to process:', (500, 100))
        #sbs = wx.StaticBoxSizer(sb, orient=wx.VERTICAL)  
        self.ChooseSignal=(wx.RadioButton(self.background, 1, label='Mvalue', pos=(500, 120), style=wx.RB_GROUP))
        self.ChooseSignal=(wx.RadioButton(self.background, 1, label='Avalue', pos=(500, 140)))
        
        wx.StaticText(self.background, 2, 'Choose the microarray design:', (500, 200))
        #sbs = wx.StaticBoxSizer(sb, orient=wx.VERTICAL)  
        self.ChooseDesign=(wx.RadioButton(self.background, 2, label='Tiled probes', pos=(500, 220), style=wx.RB_GROUP))
        self.ChooseDesign=(wx.RadioButton(self.background, 2, label='Low-density probes', pos=(500, 240)))
       
        
        
        ##img3=wx.Image( "./imgs_chipsad_gui/punto_interr_rosso.jpeg",  wx.BITMAP_TYPE_ANY).ConvertToBitmap()
        ##self.bitmap4 = wx.StaticBitmap(self.background, -1, img3, (750, 10))#, size=(400, 150))
        self.MyInpButton=wx.Button(self.background, label='Help', pos=(780, 20))
        self.MyInpButton.Bind(wx.EVT_BUTTON, self.helpButton)
        
        ##img2=wx.Image( "./imgs_chipsad_gui/punto_esclam_rosso.jpeg",  wx.BITMAP_TYPE_ANY).ConvertToBitmap()
        ##self.bitmap3 = wx.StaticBitmap(self.background, -1, img2, (20, 270))#, size=(400, 150))
        wx.StaticText(self.background, -1, 'Algorithm options:', (40, 280) )
        
        
        self.MyInpButton=wx.Button(self.background, label='Load input file', pos=(40, 310))
        self.MyInpButton.Bind(wx.EVT_BUTTON, self.inputButton)
        self.transferArea1 = wx.TextCtrl(self.background, pos=(40,340), size=wx.Size(380, 40), style = wx.TE_READONLY |  wx.TE_BESTWRAP)
        
        
        
        
        self.MyProbeButton=wx.Button(self.background, label='Load probe file' , pos=(40, 390))
        self.MyProbeButton.Bind(wx.EVT_BUTTON, self.probeButton)
        self.transferArea2 = wx.TextCtrl(self.background, style = wx.TE_READONLY  | wx.TE_BESTWRAP, pos=(40,420), size=wx.Size(380, 40))
        
        
        
        self.MyttestButton=wx.Button(self.background, label='Set the output folder' , pos=(40, 470))
        self.MyttestButton.Bind(wx.EVT_BUTTON, self.outputButton)
        self.transferArea3 = wx.TextCtrl(self.background, style = wx.TE_READONLY  | wx.TE_BESTWRAP, pos=(40,500), size=wx.Size(380, 40))
        

        wx.StaticText(self.background, -1, 'Set the sliding window width. Default 40: ', (440, 310))
        self.sc1 = wx.SpinCtrl(self.background, -1, '',  (770, 300), (60, -1))
        self.sc1.SetRange(3,200)
        self.sc1.SetValue(40)
        
        wx.StaticText(self.background, -1, 'Set the minimun probe number for CPRs. Default 3: ', (440, 380))
        self.sc2 = wx.SpinCtrl(self.background, -1, '',  (770, 370), (60, -1))
        self.sc2.SetRange(1,10)
        self.sc2.SetValue(3)
        
        wx.StaticText(self.background, -1, 'Set the computation step. Default 10: ', (440, 450))
        self.sc3 = wx.SpinCtrl(self.background, -1, '',  (770, 440), (60, -1))
        self.sc3.SetRange(1,100)
        self.sc3.SetValue(10)
        
        self.MyComputeButton=wx.Button(self.background, label='Run' , pos=(750, 500))
        self.MyComputeButton.Bind(wx.EVT_BUTTON, self.Compute)
        
        self.MyCloseButton=wx.Button(self.background, label='Close' , pos=(420, 500))
        self.MyCloseButton.Bind(wx.EVT_BUTTON, self.Close)
        
        self.timer = wx.Timer(self.background, 1)
        self.gauge = wx.Gauge(self.background, -1,  size=(250, 25), pos=(500, 550))
        
        
        self.Show()
        
       
    
    def inputButton(self, event):
        dlg = wx.FileDialog(None, "Select a file", ".",
                           "default.txt",
                           "text file (*.txt)|*.txt|all files|*.*",
                           wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            wx.MessageBox("the choosen file is " + dlg.GetFilename() +
                         "\n\npath  " + dlg.GetPath())
            
            input=dlg.GetPath()
            self.transferArea1.SetValue(input)
            input_sad=self.transferArea1.SetValue(input)
            #self.inputvar.SetValue(input)
        # Ricordatevi sempre di distruggere la finestra!
        dlg.Destroy()
        return 1
    
    def probeButton(self, event):
        dlg = wx.FileDialog(None, "Select a file", ".",
                           "default.txt",
                           "text file (*.txt)|*.txt|all files|*.*",
                           wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            wx.MessageBox("the choosen file is " + dlg.GetFilename() +
                         "\npath  " + dlg.GetPath())
            
            probe=dlg.GetPath()
            self.transferArea2.SetValue(probe)
        # Ricordatevi sempre di distruggere la finestra!
        dlg.Destroy()
        return 1
    
    def outputButton(self, event):
        dialog = wx.DirDialog(None, "Choose a directory:",style=wx.DD_DEFAULT_STYLE | wx.DD_NEW_DIR_BUTTON)
        if dialog.ShowModal() == wx.ID_OK:
            
            wx.MessageBox("The selected folder is " + dialog.GetPath() )
            output=dialog.GetPath()
            self.transferArea3.SetValue(output)
        dialog.Destroy()

        return 1
    
    def helpButton(self, event):
        
        wx.MessageBox('Help window:', 'Info', 
            wx.OK | wx.ICON_INFORMATION)
        
        return 1
    
    def Compute(self, vent):
        signal=str(self.ChooseSignal.GetValue())
        design=str(self.ChooseDesign.GetValue())
        if signal=="False":
            signal="M"
        else:
            signal="A"
        if design=="False":
            design="h"
        else:
            design="l"
        print signal
        print design
        
        width = str(self.sc1.GetValue())
        min = str(self.sc2.GetValue())
        step = str(self.sc3.GetValue())
        input_sad=str(self.transferArea1.GetValue())
        probe_sad=str(self.transferArea2.GetValue())
        output_sad=str(self.transferArea3.GetValue())
        ttestsad=" ttest_table.txt "
        comand="python ./chipSAD.py "+input_sad+" "+probe_sad+" "+ttestsad+" -o "+output_sad+ " -w "+width+" -s "+step+" -m "+min
        #+"-d "+design+" -k "+signal
        print comand
        if ctypes.sizeof(ctypes.c_voidp)==8:
            ret2=os.system(comand)
            if ret2==0:
                print "Done!"
        
    def Close(self, event):

        dial = wx.MessageDialog(None, 'Are you sure you want to quit?', 'Question',
            wx.YES_NO | wx.NO_DEFAULT | wx.ICON_QUESTION)
            
        ret = dial.ShowModal()
        
        if ret == wx.ID_YES:
            self.Destroy()
        else:
            event.Skip()

 
   
   
#define ad loop your application
MyApp= wx.App(redirect=False)
frame= MainFrame()
MyApp.MainLoop()
       
