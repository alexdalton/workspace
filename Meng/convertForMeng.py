from optparse import OptionParser

parser = OptionParser()
parser.add_option("-f", "--file", dest="filename",
                  help="list of edge files to convert for Meng's DCA", metavar="FILE")

parser.add_option("-o", "--out", dest="outCode", help="identifier for output files")

(options, args) = parser.parse_args()

edgeFiles = open(options.filename, 'r')

linkFile = open("/workspace/Meng/data/" + "link_{0}.txt".format(options.outCode), "w")
nodeFile = open("/workspace/Meng/data/" + "node_{0}.txt".format(options.outCode), "w")

genes = set()
geneNum = 1

edgeLetter = 0
letters = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']

for edgeFileName in edgeFiles:

    maxWeight = 0
    edgeFile = open(edgeFileName.strip(), "r")
    for edge in edgeFile:
        items = edge.split()
        if len(items) < 4:
            continue
        weight = float(items[2])
        if weight < 0:
            print("Error negative weight")
        if weight > maxWeight:
            maxWeight = weight
    edgeFile.close()

    edgeFile = open(edgeFileName.strip(), "r")
    for edge in edgeFile:
        items = edge.split()
        if len(items) < 4:
            continue
        g1 = items[0]
        g2 = items[1]
        weight = items[2]
        edgeType = items[3]
        genes.add(g1)
        genes.add(g2)
        linkFile.write("{0} {1} {2} {3}\n".format(g1, g2, float(weight) / maxWeight, letters[edgeLetter]))
    edgeFile.close()
    edgeLetter += 1

for gene in genes:
    nodeFile.write("{0} {1}\n".format(gene, gene[0]))

print("Length of nodeFile: {0}".format(len(genes)))

linkFile.close()
nodeFile.close()