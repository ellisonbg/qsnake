def plot1(r,y1):
    makeplot(r,[
        (y1,"b+","1"),
        ])

def plot2(r,y1,y2):
    makeplot(r,[
        (y1,"b+","1"),
        (y2,"g+","2"),
        ])

def plot3(r,y1,y2,y3):
    makeplot(r,[
        (y1,"b+","1"),
        (y2,"g+","2"),
        (y3,"kx","3"),
        ])

def makeplot(x, ys, title="Comparison", xleg=None, yleg=None, lw=2):
    from pylab import figure, plot, setp, xlabel, ylabel, legend, xticks, \
            yticks, show, semilogx, grid
    from pylab import title as tit
    f=figure()
    leg=[]
    for g in ys:
        plot(x, g[0], g[1], lw=lw)
        leg.append(g[2])
    setp(xticks()[1], fontsize=14)
    setp(yticks()[1], fontsize=14)
    semilogx()
    grid(True)
    tit(title)
    xlabel(xleg, fontsize=14)
    ylabel(yleg, fontsize=14)
    legend(leg)
    show()

def readdns(filename):
    f=open(filename)
    f.readline()
    n=[]
    for line in f:
        nr=line.split(" ")[-1]
        n.append(float(nr))
    return n

def readur(filename):
    def convert(line):
        n=[]
        line=line.split(" ")
        while line!=[]:
            while line[0]=="":
                del line[0]
            n.append(float(line[0]))
            if len(n)==500:
                break
            del line[0]
        return n
    f=open(filename)
    return convert(f.readline()),convert(f.readline()), \
        convert(f.readline()),convert(f.readline())
