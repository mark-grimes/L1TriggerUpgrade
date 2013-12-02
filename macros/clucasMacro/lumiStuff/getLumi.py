from optparse import OptionParser
from sys import exit

parser = OptionParser()

parser.add_option("-i", "--infile",
                  dest="inFileName", default=0,
                  help="input file from lumiCalc")
parser.add_option("-o", "--outfile",
                  dest="outFileName", default=0,
                  help="output file")
parser.add_option("-p", "--pixCorr",
                  action="store_true", dest="doPixCorrection", default=False,
                  help="apply the pixel luminosity correction (for 7TeV HPF)")
parser.add_option("-v", "--verbose",
                  action="store_true", dest="verbose", default=False,
                  help="turn on detailed terminal output")
parser.add_option("-n", "--nBunches",
				  dest="numBunches", default=2,
				  help="number of bunches in run")

(options, args) = parser.parse_args()

if (options.inFileName == 0) or (options.outFileName == 0):
	print "\n\t***Error: Please specify input and output files."
	print "\tFor help, run \"python getLumi.py -help\"\n"
	exit()

def pixelCorrection(rawLumi):
	#nBX=2. #CHANGE FOR DIFFERENT RUNS
	nBX=float(options.numBunches)
	print nBX
	
	val          = rawLumi
	gamma        = 1.141
	beta         = 1.
	fOrbit       = 11246.
	secondsPerLS = pow(2, 18) / fOrbit
	Ldot         = rawLumi/(nBX*secondsPerLS)
	alpha        = 0.076
	
	val     = rawLumi*gamma*beta*(1. - alpha*Ldot)
	
	alpha_1 = 0.0700
	alpha_2 = -0.0045
	val     = rawLumi*gamma*beta/(1. + alpha_1*Ldot + alpha_2*Ldot*Ldot);
	
	return val

def calcPU(intLumi):
	minBiasXS = 69.3E3 #converted to units of ub-1...
	tLumi     = 23.3570304
	fOrbit    = 11246.
	#nBX      =2.
	nBX       =float(options.numBunches)

	val = (intLumi*minBiasXS)/(tLumi*fOrbit*nBX)
	
	return val

def calcInstLumi(intLumi):
	return intLumi*pow(10,30)/23.3570304
	

f1 = open(options.inFileName, 'r')
f2 = open(options.outFileName, 'w')

for line in f1:
	lineTmp = line.split(",")
  	if options.doPixCorrection:
  	    corrLumi = pixelCorrection(float( lineTmp[6] )) #using manually applied pixel correction
   	else:
   		corrLumi = float( lineTmp[6]) #using standard luminosity correction
   	if options.verbose:
   		print "LS ", lineTmp[1]
   		print "IntLumi ", corrLumi
   		print "InstLumi", calcInstLumi(corrLumi)
   		print "PU ", calcPU(corrLumi), "\n"
	f2.write("%s\n%s\n%s\n%s\n" % (lineTmp[1], corrLumi, calcInstLumi(corrLumi), calcPU(corrLumi)))

if options.doPixCorrection:
	print "*** Ran using non-standard Pixel Corrections ***"	
print "\nCorrected Lumi and PU values printed to", options.outFileName, "\n"
