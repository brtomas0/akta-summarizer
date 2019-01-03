# TODOS: make a maping function that maps arbitrary number of points from one y-value array to 
# another array of x-values. Need this to take pH values (that are different in number from UV)
# and give every uv x-value a pH value that is linear-interpolated/approximated from a pH array

# must pip install pycorn
#      pip install matplotlib
#      pip install scipy

import os
import xml.etree.ElementTree as ET
import pycorn
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal


class Event():
    def __init__(self, element):
        # element from elementtree class
        self.EventSubType = element.get("EventSubType")
        self.EventTime = self.getText(element, "EventTime")
        self.EventTime_index = None
        self.EventVolume = self.getText(element, "EventVolume")
        self.EventVolume_index = None
        self.EventText = self.getText(element, "EventText")
        self.InstructionFeedback = self.getText(element, "InstructionFeedback")

    def getEventTimeIndex(self, time_list):
        if self.EventTime_index != None:
            return self.EventTime_index
        else:
            for i,value in enumerate(volume_data):
                if float(value) >= float(self.EventVolume):
                    self.EventTime_index = i
                    return i
            raise ValueError('Could not find %s in time_list'%value)

    def getEventVolumeIndex(self, volume_data):
        # assumes volume_data is ordered from least->greatest
        if self.EventTime_index != None:
            return self.EventTime_index
        else:
            for i,value in enumerate(volume_data):
                if float(value) >= float(self.EventVolume):
                    self.EventTime_index = i
                    return i
            # if volume is greater than data recorded, return the end of volumedata
            return(len(volume_data) - 1)

            # print(self.EventVolume, volume_data[-1])
            # raise ValueError('Could not find %s in volume data for %s'%(self.EventVolume, self.EventText))

    def getText(self, element, subject):
        # the return text from a found element will return None if none found, I want it to return a string
        s = element.find(subject).text
        if s == None:
            return ""
        else:
            return s

    def __str__(self):
        return(self.EventVolume.ljust(12) + self.EventText)

def XMLTesting():
    pycorn_file = pycorn.pc_uni6("VenB145-028\\1.0(636803223998414676).zip")
    pycorn_file.load()

    chrom_1_tree = ET.fromstring(pycorn_file["Chrom.1.Xml"])

    load_volume = 0
    peak_volume = 0
    in_sample_app_phase = False
    spacer = ""

    for child in chrom_1_tree.iterfind("./EventCurves//Event"): 
        print_message = False

        event = Event(child)

        if event.EventSubType == "BlockStart":
            spacer +=  "\t|"
            print_message = True

            if "Sample Application" in event.EventText:
                in_sample_app_phase = True
            elif in_sample_app_phase and "Start frac" in event.EventText:
                load_volume -= float(child.find("EventVolume").text)
            elif in_sample_app_phase and "Stop frac" in event.EventText:
                load_volume += float(child.find("EventVolume").text)
            
        elif event.EventSubType == "BlockEnd":
            spacer = spacer[:-2]
            print_message = False

            if in_sample_app_phase and len(spacer)/2<2:
                in_sample_app_phase = False
            
        elif "Peak start" in event.EventText:
            peak_volume -= float(child.find("EventVolume").text)
            print_message = True
        elif "Peak end" in event.EventText:
            peak_volume += float(child.find("EventVolume").text)
            print_message = True

        elif "Frac" in event.EventText:
            print_message = False 

        if print_message:
            message = child.find("EventText").text
            volume = float(child.find("EventVolume").text)
            print(str(volume).ljust(12) + spacer + message[:message.find("(Issued)")])

    print("Load Volume: %0.2f"%load_volume)
    print("Peak Volume: %0.3f"%peak_volume)

def getEventIndices(chrom_1_tree, volume_data):
    num_volume_data = len(volume_data) - 1
    indices = { 'sampleApp': [0, num_volume_data],
                'FT': [0, num_volume_data],
                'elution': [0, num_volume_data],
                'peak': [0, num_volume_data],
                'cleaning': [0, num_volume_data]}

    depth = 0
    save_depth = [float("inf")]
    save_key = [None]
    for child in chrom_1_tree.iterfind("./EventCurves//Event"):
        event = Event(child)
        # print("\t"*depth + event.EventText)

        if event.EventSubType == "BlockStart":
            depth += 1

            if "Phase Elution" in event.EventText:
                indices['elution'][0] = event.getEventVolumeIndex(volume_data)
                save_depth.append(depth)
                save_key.append('elution')

            elif "Sample Application" in event.EventText or "Sample Load" in event.EventText:
                indices['sampleApp'][0] = event.getEventVolumeIndex(volume_data)
                save_depth.append(depth)
                save_key.append('sampleApp')

            elif save_key[-1] == "sampleApp":
                if "Start frac" in event.EventText:
                    indices['FT'][0] = event.getEventVolumeIndex(volume_data)
                if "Stop frac" in event.EventText:
                    indices['FT'][1] = event.getEventVolumeIndex(volume_data)

        elif event.EventSubType == "BlockEnd":
            if save_depth[-1] == depth and save_key[-1] == 'elution':
                indices['elution'][1] = event.getEventVolumeIndex(volume_data)
                save_depth.pop()
                save_key.pop()
            elif save_depth[-1] == depth and save_key[-1] == 'sampleApp':
                indices['sampleApp'][1] = event.getEventVolumeIndex(volume_data)
                save_depth.pop()
                save_key.pop()
            depth -= 1

        elif "NaOH Clean" in event.EventText:
            indices['cleaning'][0] = event.getEventVolumeIndex(volume_data)
        elif "Neutralization" in event.EventText:
            indices['cleaning'][1] = event.getEventVolumeIndex(volume_data)

        elif "Peak start" in event.EventText:
            indices['peak'][0] = event.getEventVolumeIndex(volume_data)
        elif "Peak end" in event.EventText:
            indices['peak'][1] = event.getEventVolumeIndex(volume_data)
        elif indices['peak'][1] == num_volume_data and "Stop peak" in event.EventText:
            # added in case the peak goes to end of fractionation
            indices['peak'][1] = event.getEventVolumeIndex(volume_data)

    return indices

def getData(filepath):
    #print(file)
    pycorn_file = pycorn.pc_uni6(filepath)
    pycorn_file.load()
    chrom_1_tree = ET.fromstring(pycorn_file["Chrom.1.Xml"])

    volume_data = pycorn_file["Chrom.1_1_True"]["CoordinateData.Volumes"]
    UV280_data = pycorn_file["Chrom.1_1_True"]["CoordinateData.Amplitudes"]
    pH_data = pycorn_file["Chrom.1_7_True"]["CoordinateData.Amplitudes"] #is taken at different time points

    indices = getEventIndices(chrom_1_tree, volume_data)

    return chrom_1_tree, volume_data, UV280_data, indices

def alignmentIndex(location, volume_data, UV280_data, indices):
    loc_split = location.split(" ")

    if len(loc_split) == 1:
        if location in indices:
            return indices[loc_split[0]][0]
        else:
            print("first key must be one of these: {%s}"%indices.keys())
            return 0

    if loc_split[0] not in indices:
        print("first key must be one of these: {%s}"%indices.keys())
        return 0
    elif loc_split[1] not in ('peak', 'begin', 'end'):
        print("second key must be one of these: {%s}"%list(('peak', 'begin', 'end')))
        return indices[loc_split[0]][0]
    elif loc_split[1] =='peak':
        return UV280_data.index(max(UV280_data[indices[loc_split[0]][0]:indices[loc_split[0]][1]]))
    elif loc_split[1] == 'begin':
        return indices[loc_split[0]][0]
    elif loc_split[1] == 'end':
        return indices[loc_split[0]][1]
    else:
        return indices[loc_split[0]][0]

def smoothSpikes(data, window = None):
    if window == None:
        window = int(0.01 * len(data))
    pass

def analyzeGEFiles(input_folder, output_folder, volume_align, UV_align):
    if input_folder not in os.listdir():
        print('Input Folder "%s" not in current directory'%input_folder)
        return

    file_list = os.listdir(input_folder)
    number_of_files = len(file_list)

    for filenumber, file in enumerate(file_list):
        try:
            sample_name = file[:file.index("(")]
        except:
            sample_name = file[:file.index(".")]

        with open(output_folder + "\\" + file[:-4] + "_output.csv", "w") as ofile:
            chrom_1_tree, volume_data, UV280_data, indices = getData(input_folder + "\\" + file)
            print(indices)

            #Shifting to align the peaks and a value in the elution after it has stabilized
            FT_volume = volume_data[indices['sampleApp'][1]] - volume_data[indices['sampleApp'][0]]
            peak_volume = volume_data[indices['peak'][1]] - volume_data[indices['peak'][0]]

            volume_align_index = alignmentIndex(volume_align, volume_data, UV280_data, indices)
            UV280_align_index = alignmentIndex(UV_align, volume_data, UV280_data, indices)

            print("%s,%.3f,%.3f"%(sample_name, FT_volume, peak_volume))
            
            # Write the file header
            ofile.write("Filename:,%s\n"%file)
            ofile.write("Elution Phase Begin (mL):,%.3f\n"%volume_data[indices['elution'][0]])
            ofile.write("Load/FT Volume (mL):,%.3f\n"%FT_volume)
            ofile.write("Elution Peak Volume (mL):,%.3f\n"%peak_volume)

            # write the volume,UV data to the file
            for i in range(len(UV280_data)):
                ofile.write("%s,%s\n"%(volume_data[i] - volume_data[indices['peak'][0]], UV280_data[i]))

            # Make the figure
            xdata = [x - volume_data[volume_align_index] for x in volume_data]
            ydata = [y - UV280_data[UV280_align_index] for y in UV280_data]

            ydata_smoothed = signal.medfilt(ydata, int(len(ydata)/300) * 2 + 1)
            ydata_smoothed = signal.savgol_filter(ydata_smoothed, int(len(ydata_smoothed)/400) * 2 + 1, 2)

            line_color = "#%02x90%02x"%(int((number_of_files-filenumber-1)*255/(number_of_files-1)), int((filenumber)*255/(number_of_files-1)))

            plt.plot(xdata, ydata, label = sample_name, color = line_color)
            plt.plot(xdata, ydata_smoothed, color = line_color, linewidth = 3)

    # Set graphing variables for 1920/1080 viewing
    plt.xlabel("Volume (L)", fontsize = 14)
    plt.ylabel("UV 280 Absorbance (Au)", fontsize = 14)
    plt.legend(loc = "upper left", fontsize = 18, markerscale = 10)
    plt.title(input_folder + " UV-280 Elution Peaks", fontsize = 24)
    plt.subplots_adjust(left = 0.04, right = 0.99, bottom = 0.05, top = 0.95)
    plt.show()

def main():
    input_folder = "input_folder"
    output_folder = "output_folder"
    volume_align = "elution end"
    UV_align = "elution end"

    analyzeGEFiles(input_folder, output_folder, volume_align, UV_align)

def openFunction():
    for i in open("output.csv").readlines()[0:5]: print(i.strip())

def testing():
    a = 3
    b = 255
    for i in range(a):
        print(i*b/(a-1))

if __name__ == '__main__':
    main()
    