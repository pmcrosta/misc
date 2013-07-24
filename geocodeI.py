###########################
## This program defines functions that utilize the google and yahoo geocoder APIs to
## convert addresses into coordinates
## Written by: PeterM Crosta
## 2/16/2009
#
#INPUT FILE MUST BE CSV ADDRESSES: "STREET, CITY, STATE, ZIPCODE"
## WARNING: THIS CODE PROBABLY DOES NOT WORK ANYMORE GIVEN THAT THE APIs HAVE CHANGED.
## POSTING HERE JUST TO KEEP TRACK OF VARIOUS PROGRAMS WRITTEN OVER THe YEARS
##########################

#import various modules
import sys, string
import urllib
from urllib import *
from xml.dom import minidom
import xml
import time

def google(fname, oname, zipo=0):
    """Use google to geocode. Need input and output csv files.
    Input must be CSV address. zipo=0 means dont run in zip code only mode"""

    #Open csv file for reading and one for writing
    f = open(fname)
    o = open(oname, 'w')

    #This is my API. You should get your own.
    api = ""

    #Loop over the lines in the input file
    for line in f:
        #rest a bit so geocoder is not flooded with requests
        time.sleep(0.5)

        #make some cosmetic changes to the line that is read in
        lines = line.replace('\t', '+')
        lines = lines.replace('\n', '')
        lines = lines.replace(' ', '+')
        lines = lines.replace('#', '+')

        #If running in zipcode only mode, just use the zipcode as input. Otherwise, use entire address
        if zipo == 1:
            code = lines.split(',')[3]
            #Add leading zero to zips that begin with 0"
            if len(code) == 4:
                code = "0"+code
            #Get xml
            site = urllib.urlopen("http://maps.google.com/maps/geo?q="+code+"&output=xml&key="+api)
        else:
            #Get xml
            site = urllib.urlopen("http://maps.google.com/maps/geo?q="+lines+"&output=xml&key="+api)

        #This creates a string of the xml
        x = site.read()

        #assume only one response unless otherwise
        multi = 1
        try:
            #create xml document object; very error prone
            kml = xml.dom.minidom.parseString(x)

            #If more than one response, count how many
            if len(kml.getElementsByTagName('AddressDetails')) > 1:
                multi = len(kml.getElementsByTagName('AddressDetails'))

            #pull out accuracy information.
            #more info http://code.google.com/apis/maps/documentation/reference.html#GGeoAddressAccuracy
            acc = kml.getElementsByTagName('AddressDetails')[0].attributes["Accuracy"].value

            #pull out lat and long
            (lng, lat) = kml.getElementsByTagName("Point")[0].getElementsByTagName("coordinates")[0].firstChild.nodeValue.split(",")[0:2]

            #extract rest of matched address depending on mode and other factors        
            if zipo==0:
                (place, city, statezip, cntry) = kml.getElementsByTagName("Placemark")[0].getElementsByTagName("address")[0].firstChild.nodeValue.split(',')[0:4]
            else:
                place = " "
                if len(kml.getElementsByTagName("Placemark")[0].getElementsByTagName("address")[0].firstChild.nodeValue.split(',')) == 3:
                    (city, statezip, cntry) = kml.getElementsByTagName("Placemark")[0].getElementsByTagName("address")[0].firstChild.nodeValue.split(',')[0:3]
                else:
                    (statezip, cntry) = kml.getElementsByTagName("Placemark")[0].getElementsByTagName("address")[0].firstChild.nodeValue.split(',')[0:2]
                    city = " "
            #cleaning
            state = statezip.strip().split(' ')[0]
            zippy = statezip.strip().split(' ')[1]

            v=acc, ",", str(multi), ",", lat, ",", lng, ",", place, ",", city, ",", state, ",", zippy, "\n"

            #write output to file
            o.writelines(v)

        except:
            #or write error to file
            o.write("Error\n")

    #close files        
    f.close()
    o.close()

def yahoo(fname, oname):
    """Use Yahoo to geocode. Need input and outfiles as arguments.
    Input file must be CSV in order: street, city, state, zip"""

    #Open csv file for reading and one for writing
    f = open(fname)
    o = open(oname, 'w')

    #This is my API. You should get your own.
    api = ""

    #Loop over the lines in the input file
    for line in f:
        #Rest a little so geocoder isnt overloaded
        time.sleep(0.5)
        
        #make some cosmetic changes to the line that is read in
        lines = line.replace('\t', '+')
        lines = lines.replace('\n', '')

        nocomma = lines.split(',')
        street, city, state, zipcode = nocomma[0].strip(), nocomma[1].strip(), nocomma[2].strip(), nocomma[3].strip()  
        street = street.replace(' ', '+')
        city = city.replace(' ', '+')

        try:        
            #Get xml
            site = urllib.urlopen("http://local.yahooapis.com/MapsService/V1/geocode?appid="+api+"&street="+street+"&city="+city+"&state="+state+"&zip="+zipcode)
            #Create string of xml
            x = site.read()
            
            try:
                #create xml document object; very error prone
                kml = xml.dom.minidom.parseString(x)

                #assume one response unless otherwise
                multi = 1
                
                if len(kml.getElementsByTagName('Error')) == 0:
                    if len(kml.getElementsByTagName('Result')) > 1:
                        multi = len(kml.getElementsByTagName('Result'))

                    #Extract address pieces from xml
                    acc = kml.getElementsByTagName('Result')[0].attributes["precision"].value
                    lat = kml.getElementsByTagName("Result")[0].getElementsByTagName("Latitude")[0].firstChild.nodeValue
                    lng = kml.getElementsByTagName("Result")[0].getElementsByTagName("Longitude")[0].firstChild.nodeValue
                    place = kml.getElementsByTagName("Result")[0].getElementsByTagName("Address")[0].firstChild.nodeValue
                    cityy = kml.getElementsByTagName("Result")[0].getElementsByTagName("City")[0].firstChild.nodeValue
                    statey = kml.getElementsByTagName("Result")[0].getElementsByTagName("State")[0].firstChild.nodeValue
                    zipy = kml.getElementsByTagName("Result")[0].getElementsByTagName("Zip")[0].firstChild.nodeValue
                
                    v = acc, ",", str(multi), ",", lat, ",", lng, ",", place, ",", cityy, ",", statey, ",", zipy, "\n"

                    #write output to file
                    o.writelines(v)

                else:
                    #or write error to file
                    o.write("Error\n")
                    
            except:
                #or write error to file
                o.write("Error\n")
                
        except:
            #or write error to file
            o.write("Error\n")

    #close files
    f.close()
    o.close()
