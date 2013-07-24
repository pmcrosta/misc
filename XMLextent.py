###########################
## This program extracts four coordinates from an XML file and creates
## four coordinate pairs.
## Written by: Peter M Crosta
## 8/28/2009
## For use with ArcGIS GeoProcessor
##########################

#import various modules
import sys, string, os, arcgisscripting
from xml.dom import minidom
import xml
import time
from Tkinter import Tk
import tkFileDialog

#set up the geoprocessor
gp = arcgisscripting.create(9.3)
gp.SetProduct("ArcInfo")
gp.toolbox = "management"
#gp.Workspace = sys.path[0]
gp.Overwriteoutput = 1

#selecting output location
toplevel=Tk()
toplevel.withdraw()
sname = tkFileDialog.askdirectory(title='Select directory for output shape files')
toplevel.destroy()

def xmlit():
    """no inputs"""
    
    #Open xml file for reading
    #fname = "CommercialMapOfChina.xml"
    toplevel=Tk()
    toplevel.withdraw()
    fname=tkFileDialog.askopenfilename(title='Please select an XML', filetypes=[('XML files','*.xml')])
    toplevel.destroy()
    f = open(fname)

    lines = f.read()

    try:
        #create xml document object; very error prone
        kml = xml.dom.minidom.parseString(lines)

        #extract xml data
        westbc = kml.getElementsByTagName("bounding")[0].getElementsByTagName("westbc")[0].firstChild.nodeValue[1:]
        eastbc = kml.getElementsByTagName("bounding")[0].getElementsByTagName("eastbc")[0].firstChild.nodeValue[1:]
        northbc = kml.getElementsByTagName("bounding")[0].getElementsByTagName("northbc")[0].firstChild.nodeValue[1:]
        southbc = kml.getElementsByTagName("bounding")[0].getElementsByTagName("southbc")[0].firstChild.nodeValue[1:]

        title = kml.getElementsByTagName("title")[0].firstChild.nodeValue
        link = kml.getElementsByTagName("onlink")[0].firstChild.nodeValue
        thumb = kml.getElementsByTagName("browsen")[0].firstChild.nodeValue
        
        #create extent coordinates
        ws = (westbc, southbc)
        wn = (westbc, northbc)
        en = (eastbc, northbc)
        es = (eastbc, southbc)

        print "Title-Link-Thumbnail"
        print title+"\n"+link+"\n"+thumb+"\n"
        print "coordinates\n"
        print ws, wn, en, es
        
    except:
        print "something went wrong."


    #close files        
    f.close()

    
    ##create new feature class
    gp.CreateFeatureClass(sname, os.path.basename(fname)[:-3]+"shp", "POINT")
    fc = sname+"\\"+os.path.basename(fname)[:-3]+"shp"

    #create point opbject
    pnt = gp.CreateObject("Point")
    #instantiate cursor in new feature class
    cur = gp.InsertCursor(fc)
    #create describe object
    desc = gp.Describe(fc)
    #pull out field that describes shape
    shpFld = desc.ShapeFieldName

    #add four points to new feature class
    pnt.X, pnt.Y = ws
    feat = cur.NewRow()
    feat.SetValue(shpFld, pnt)
    cur.InsertRow(feat)
    pnt.X, pnt.Y = wn
    feat = cur.NewRow()
    feat.SetValue(shpFld, pnt)
    cur.InsertRow(feat)
    pnt.X, pnt.Y = en
    feat = cur.NewRow()
    feat.SetValue(shpFld, pnt)
    cur.InsertRow(feat)
    pnt.X, pnt.Y = es
    feat = cur.NewRow()
    feat.SetValue(shpFld, pnt)
    cur.InsertRow(feat)

    del pnt, cur, feat, desc, shpFld, fc, ws, wn, en, es, title, link, thumb
