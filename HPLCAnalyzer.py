#!/bin/python

from optparse import OptionParser
from scipy import integrate, signal
import numpy
import matplotlib.pyplot as plt
from matplotlib.text import OffsetFrom
import matplotlib

class peak:
    def __init__(self, chrom, startTime, endTime, color='b' ):
        self.ritention_time = None
        self.area = None
        self.chrom = chrom
        self.integration_start = self.time_to_index(startTime)
        self.integration_end = self.time_to_index(endTime)
        self.startP = [self.chrom.chrom[0][self.integration_start],
                       self.chrom.chrom[1][self.integration_start]]
        self.endP = [self.chrom.chrom[0][self.integration_end],
                     self.chrom.chrom[1][self.integration_end]]
        self.color = color
        peakChrom = self.chrom.chrom[1][self.integration_start:self.integration_end]
        self.time = self.chrom.chrom[0][peakChrom.argmax(0) + self.integration_start]
        self.peakMaximum = self.chrom.chrom[1][peakChrom.argmax(0) + self.integration_start]

    def time_to_index(self, time):
        timedifference = -1
        time_found = 0
        for t in self.chrom.chrom[0]:
            if abs(t - time) < timedifference or timedifference < 0:
                timedifference = abs(t - time)
                time_found = t

        return self.chrom.chrom[0].index(time_found)

    def baseline( self ):
        return ( self.endP[1] + self.startP[1] ) * ( self.endP[0] - self.startP[0] ) / 2

    def calculateArea( self ):
        peak = []
        peak.append([])
        peak.append([])
        peak[0] = self.chrom.chrom[0][self.integration_start:self.integration_end]
        peak[1] = self.chrom.chrom[1][self.integration_start:self.integration_end]
        baseline_correction = self.baseline()
        self.chrom.MPLax.fill(peak[0], peak[1], self.color, alpha=0.3)
        self.chrom.MPLfig.canvas.draw()
        self.area = integrate.trapz( peak[1], peak[0] , 0.01 ) - baseline_correction

class chromatogram:
    def __init__(self, file = None, normalized = False ):
        self.name = ""
        self.chrom = []
        self.chrom.append([])
        self.chrom.append([])

        self.P = []
        self.peak_annotations = []
        self.__COLORS__ = ['r', 'b', 'y', 'g']
        self.deadtime = 0.0

        if file is not None:
            self.name = file.split('.')[0]
            self.importChromatogramFromXY( file )
            if normalized == True:
                self.normalizeChromatogram()
            self.plotChromatogram()

    def plotChromatogram( self ):
        self.font = {'family' : 'serif', 'size': 25}
        matplotlib.rc( 'font', **self.font )
        self.MPLfig = plt.figure(1)
        self.MPLax = self.MPLfig.add_subplot(111)
        self.MPLax.set_xlabel('time (min)')
        self.MPLax.get_yaxis().set_visible(False)
        self.MPLax.spines["top"].set_visible(False)
        self.MPLax.spines["left"].set_visible(False)
        self.MPLax.spines["right"].set_visible(False)
        self.MPLax.set_xlim([0, self.chrom[0][len(self.chrom[0])-1]])
        self.MPLax.set_ylim([0, 1.15 * numpy.amax( self.chrom[1] )])

        self.clickevent = self.MPLfig.canvas.mpl_connect('button_press_event', self.__onclick__add_peak )
        self.Revent = self.MPLfig.canvas.mpl_connect('key_press_event', self.__on_R_key__)
        self.Pevent = self.MPLfig.canvas.mpl_connect('key_press_event', self.__on_A_key__)

        self.Click1 = None
        self.MPLax.plot(self.chrom[0], self.chrom[1], color = 'k', lw=0.5)

    def importChromatogramFromXY( self, fn ):
        infile = open(fn)
        data = infile.read()
        infile.close()
        for row in data.split('\n'):
            if row != "":
                coord = row.split(',')
                self.chrom[0].append(float(coord[0]))
                self.chrom[1].append(float(coord[1]))

    def normalizeChromatogram( self ):
        self.chrom = [self.chrom[0],
                      numpy.multiply( self.chrom[1], 1/numpy.amax( self.chrom[1] ) )]
        self.chrom = [self.chrom[0],
                      numpy.subtract( self.chrom[1], numpy.amin( self.chrom[1] ) )]

    def addPeak( self, start, end, color=None ):
        if color is None:
            color = self.__COLORS__[(len(self.P)) % len(self.__COLORS__)]

        self.P.append( peak( self, start, end, color ) )
        self.P[len(self.P)-1].calculateArea()
        self.refreshLabels()

    def refreshLabels( self ):
        for i in self.peak_annotations:
            i.remove()
            self.peak_annotations.remove(i)

        #Generate label for t0 = deadtime
        if self.deadtime != 0:
            deadtime_index = 0
            for i in range( len( self.chrom[0] ) ):
                if abs( self.deadtime - self.chrom[0][i] ) < abs( self.deadtime - self.chrom[0][deadtime_index] ):
                    deadtime_index = i
            self.peak_annotations.append( self.MPLax.annotate("t$_0$ = {:.1f}".format(self.deadtime),
                                xy=(self.deadtime, self.chrom[1][deadtime_index]), xycoords="data",
                              xytext=(self.deadtime, 0.2 * numpy.amax( self.chrom[1] ) + self.chrom[1][deadtime_index] ),
                              # xytext is offset points from "xy=(0.5, 0), xycoords=an1"
                              va="top", ha="center",
                              arrowprops=dict(arrowstyle="]-", lw=0.8, facecolor='black')))

        for i in range( len(self.P) ):
            self.generateLabels( i )

    def generateLabels( self, peak_index ):
        #Area label
        totalarea = 0.0
        for p in self.P:
            totalarea += p.area
        areaText = "{:.1f}%".format(self.P[peak_index].area / totalarea * 100 )
        self.peak_annotations.append( self.MPLax.annotate(areaText, xy=(self.P[peak_index].time, 0.8 * self.P[peak_index].peakMaximum), xycoords="data",
                          xytext=(self.P[peak_index].time + 4, self.P[peak_index].peakMaximum),
                          # xytext is offset points from "xy=(0.5, 0), xycoords=an1"
                          va="top", ha="center",
                          bbox=dict(boxstyle="round4", fc="w", lw=0.8),
                          arrowprops=dict(arrowstyle="-|>", lw=0.8, connectionstyle="angle3,angleA=90",facecolor='black')))

        self.peak_annotations.append(self.MPLax.text(self.P[peak_index].time, 0.05 + self.P[peak_index].peakMaximum,
                                                    "{:.1f}".format( self.P[peak_index].time - self.deadtime),
                                                    va="center", ha="center"))
        self.MPLfig.canvas.draw()


    def __onclick__add_peak( self, click ):
        if self.Click1 is None:
            self.Click1 = click.xdata
        else:
            if click.xdata > self.Click1:
                self.addPeak( self.Click1, click.xdata )
            else:
                self.addPeak( click.xdata, self.Click1 )


            self.Click1 = None

    def __onclick__dead_time( self, click ):
        self.deadtime = click.xdata
        self.refreshLabels()
        click.canvas.mpl_disconnect( self.clickevent )
        self.clickevent = click.canvas.mpl_connect('button_press_event', self.__onclick__add_peak )

    def __on_R_key__( self, event ):
        if event.key == 'r':
            event.canvas.mpl_disconnect( self.clickevent )
            self.clickevent = event.canvas.mpl_connect('button_press_event', self.__onclick__dead_time )

    def printed_report( self ):
        print("**************************")
        print("Dead time: {:.2f} min".format(self.deadtime) )
        if len(self.P) == 2:
            print("e.e. = {:.2f}".format(((self.P[0].area-self.P[1].area)/(self.P[0].area+self.P[1].area)*100)))

        totalarea = 0.0
        for p in self.P:
            totalarea += p.area
        for i in self.P:
            print("Peak {:d}".format( self.P.index(i) ) )
            print("\tArea = {:.2f} ({:.1f}%)".format(i.area, i.area/totalarea*100))
            print("\tRetention time = {:.2f}".format(i.time - self.deadtime) )
        print("**************************")

    def saveimage( self ):
        self.MPLfig.savefig("{}.svg".format(self.name))

    def __on_A_key__( self, event ):
        if event.key == 'a':
            self.printed_report()
            self.saveimage()


parser = OptionParser()
parser.add_option("-f", "--file", dest="filename",
                  help="open this chromatogram", metavar="FILE")

(options, args) = parser.parse_args()

chrom = chromatogram(file=options.filename, normalized=True)
plt.show()
