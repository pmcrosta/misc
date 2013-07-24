##Script Name: Change String Fields into Float Fields for selected shapefiles and fields for CDC China data
##Created By: Peter M Crosta
##Date: 31/10/2008

'''
##This script takes CDC China township datasets as inputs and converts the string fields
##that should be numeric to float. It also creates a flag variable that is set equal to 1 for each
##observation that was missing but converted to a float of zero.

New Rules: This program hangs on the last merge command if run as an ArcGIS script. It must
be run in Python instead. You can only do the conversion for files in one folder at a time.
'''

#manage inputs
#This script takes parameters for input datasets (multivalue input) and an output workspace where the final products will reside
#inputs = gp.GetParameterAsText(0)
#datasets = string.split(inputs,";")
#outputws = gp.GetParameterAsText(1)

#the following used for building and testing; do not edit
#d = "D:/MyDocs/Yunnan/townshipl2.shp"

######################################################################
###########################PROGRAM  INPUTS############################
##################Modify this block of code, nothing else#############
######################################################################

#EXAMPLE
#inputws = "D:/MyDocs/Yunnan"
#outputws = "D:/MyDocs/Yunnan/NewTwp"
#inputs = "townshipl2.shp;townshipl1.shp;townshipa.shp"
#
#Enter output folder here. It must be created already.:
outputws = "D:/MyDocs/jt2118/Gonsu/NewTwp"
#Enter input folder path here:
inputws = "D:/MyDocs/jt2118/Gonsu"
#Enter list of shapefiles to convert here. separate file names by semi-colon and no spaces:
inputs = "townshipa.shp"
######################################################################
######################################################################

#Import modules
import arcgisscripting, sys, os, string

#from time import *
#print strftime("Start: %d %b %Y %H:%M:%S", localtime())

#create geoprocessor object
gp = arcgisscripting.create()

#Split input datasets in tuple
datasets = string.split(inputs,";")

#load data management toolbox
gp.toolbox = "management"
print "Preprocessing complete"

#if d:
for d in datasets:
    print "Processing... "+d
    #Create output featureclass with schema of original file
    desc=gp.describe(inputws+"/"+d)
    gp.CreateFeatureClass(outputws, desc.name, desc.shapetype, inputws+"/"+d)
    gp.CreateFeatureClass(outputws, 'f_'+desc.name, desc.shapetype)
    
    #Create copy of original file in output WS so original data is not modified
    gp.copy_management(inputws+"/"+d, outputws + "/t_" + desc.name)
    dataset = outputws + "/t_" + desc.name
    newset = outputws + "/f_" + desc.name
    
    #Get list of string fields and then iterate until we are at the first one to convert to float
    fields = gp.ListFields(inputws+"/"+d, "", "string")
    field = fields.next()

    while field.name <> 'ENAME':
        field = fields.next()
    field = fields.next()

    print "Enter Field Loop"
    #a, b, and c are used to index lists of fields that need to be converted (fields is already limited to strings)
    a,b,c,h =[], [], [], 1
    while field:
        #Create list containing string field names to convert
        a.append(field.name)
        #Create new temporary field name and flag field name
        b.append(field.name+'_t')
        c.append(field.name+'_f')
        #Add temporary field to dataset, perm flag field, and flag field to output dataset
        gp.addfield(dataset, field.name+'_t', "float")
        gp.addfield(newset, field.name+'_f', "short")
        gp.addfield(outputws+"/"+desc.name, field.name+'_f', "short")
        #Create string list of fields to include in cursor call
        if h == 1:
            flist = field.name+"; "+field.name+'_t'
            h = 0
        else:
            flist = flist+"; "+field.name+"; "+field.name+'_t'
        field=fields.next()

    print "Exit Field Loop"
    #Define an update cursor that focuses on the two fields that will be used for conversion
    #also define insert cursor for FLAG dataset
    cursor = gp.UpdateCursor(dataset, "", "", flist)
    newRow=cursor.next()
    cursor2 = gp.InsertCursor(newset)

    print "Cursors Instantiated "    
    #loop over the rows and convert each string to float as long as it is not missing. If it is missing, the flag indicator is set to 1
    while newRow:
    #loop over fields
        for i in range(0,len(a)):         
            x = newRow.GetValue(a[i])
            newRow2=cursor2.NewRow()
            if x.isdigit():
                newRow.SetValue(b[i], float(x))
                newRow2.SetValue(c[i], 0)
            else:
                newRow2.SetValue(c[i], 1)

        #Update and select next row
        cursor.UpdateRow(newRow)
        newRow=cursor.next()
        cursor2.InsertRow(newRow2)

    #delete cursor to remove lock on dataset
    del cursor, newRow, x, cursor2, newRow2
    print "Cursors Deleted"
    for j in range(0,len(a)):
        #delete original string field
        gp.deletefield(dataset, a[j])
        #add a float field of the same name (that was just deleted) right back to the end of dataset
        gp.addfield(dataset, a[j], "float")                
        #Make copy of recently converted temp field to new float field
        oldfield = "!"+b[j]+"!"
        gp.CalculateField_management(dataset, a[j], oldfield, "PYTHON")
        #delete temporary field
        gp.deletefield(dataset, b[j])

    #Join flag fields to rest of data.
    gp.joinfield_management(dataset, "FID", newset, "FID")
    gp.deletefield(dataset, 'Id')
    
    #Create fieldmappings object from newly created dataset that has original field order plus flag fields at the end
    #This block of fieldmappings code serves to protect the field order of the original shapefile
    fieldmappings = gp.CreateObject("FieldMappings")
    fieldmappings.AddTable(outputws+"/"+desc.name)

    #typer is a list of fields in fieldmappings; iterate through to the first field that can be modified
    typer = fieldmappings.fields
    changer=typer.next()

    #This loop goes through the fields in the final output table and created a fieldmap for each one so that data are merged correctly
    #First a Fieldmap object is created out of thin air. Then the object gets the input field properties from the corresponding field in the dataset
    #This is the strangest step in this script, and I don't totally understand how it works. But it does.
    while changer:
        fieldmap = gp.CreateObject("FieldMap")
        fieldmap.addinputfield(dataset, changer.name)
        fieldmappings.replacefieldmap(fieldmappings.findfieldmapindex(changer.name), fieldmap)
        del fieldmap
        changer=typer.next()

    del typer, changer
    print "Fieldmaps done"
    
    #now merge new data into a dataset that has the original field order and use fieldmappings to point everything in the right place
    #Also, delete dataset that fieldmappings came from first and delete temporary dataset.
    gp.delete_management(outputws+"/"+desc.name)
    gp.merge(dataset, outputws+"/"+desc.name, fieldmappings)
    gp.delete_management(outputws+"/t_"+desc.name)
    gp.delete_management(outputws+"/f_"+desc.name)
    print outputws+"/"+desc.name+" has been created"

    del a, b, c, h, i, j, field, oldfield, desc, flist, fieldmappings, fields
    
del gp, inputs, outputws, datasets, dataset
