import numpy as np
import matplotlib.pyplot as plt
import sys

def getColor(id):
    cmap = plt.get_cmap('tab20').colors
    return cmap[id%len(cmap)]

def darker(color):
    import matplotlib.colors as mc
    import colorsys
    c = colorsys.rgb_to_hls(*mc.to_rgb(color))
    return colorsys.hls_to_rgb(c[0],0.75*c[1],c[2])

def minmax(data):
    minD = np.min(data)
    maxD = np.max(data)
    rangeD = maxD - minD
    if (abs(minD) < 0.1*rangeD):
        minD = -0.1*rangeD
    elif (minD < 0):
        minD = 1.1*minD
    else:
        minD = 0.9*minD
    if (abs(maxD) < 0.1*rangeD):
        maxD = 0.1*rangeD
    elif (maxD < 0):
        maxD = 0.9*maxD
    else:
        maxD = 1.1*maxD
    return minD, maxD

def genfromcsv(prefix):
    dataV = np.loadtxt("vert.csv",delimiter=",")
    dataE = np.loadtxt("edge.csv",delimiter=",")
    dataC = np.loadtxt("cell.csv",delimiter=",")
    dataEV = {}
    dataCE = {}
    with open("edgeVertices.csv", "r") as fh:
        for line in fh:
            line = line.strip()
            parts = line.split(", ")
            eId = int(parts[0])
            nbV = int((len(parts)-1)/2)
            dataEV[eId] = {"vertId": np.array(parts[1:nbV+1],dtype=np.double), "vertOrientation": np.array(parts[nbV+1:],dtype=np.double)}
    with open("cellEdges.csv", "r") as fh:
        for line in fh:
            line = line.strip()
            parts = line.split(", ")
            eId = int(parts[0])
            nbV = int((len(parts)-2)/2)
            dataCE[eId] = {"flat": bool(parts[1]), "edgeId": np.array(parts[2:nbV+2],dtype=int), "edgeOrientation": np.array(parts[nbV+2:],dtype=np.double)}
    for mId in np.unique(dataV[:,0]):
        locDataV = dataV[dataV[:,0]==mId,1:]
        locDataE = dataE[dataE[:,0]==mId,1:]
        locDataC = dataC[dataC[:,0]==mId,1:]
        # Plot with ids
        fig, ax = plt.subplots()
        fig.suptitle("Map Id {}".format(int(mId)))
        fig.set_size_inches(10,10)
        ax.set_xlim(*minmax(locDataV[:,0]))
        ax.set_ylim(*minmax(locDataV[:,1]))
        ax.set_aspect('equal')
        ax.set_title('Connectivity')
        for vert in locDataV[:]:
            ax.text(vert[0],vert[1],int(vert[2]),color="black")
        for edge in locDataE[:]:
            ax.text(edge[0],edge[1],int(edge[2]),color="coral")
        for cell in locDataC[:]:
            ax.text(cell[0],cell[1],int(cell[2]),color="lime")
        fig.savefig(prefix+"Connectivity{}.png".format(int(mId)))
        # Plot with orientation
        ax.clear()
        ax.set_title('Orientation')
        for cell in locDataC[:]:
            cId = int(cell[2])
            eIds = list(dataCE[cId]["edgeId"])
            vIds = [dataEV[eIds[0]]["vertId"][0],dataEV[eIds[0]]["vertId"][1]]
            eIds.pop(0)
            while( vIds[0] != vIds[-1]):
                for iE in eIds:
                    if (dataEV[iE]["vertId"][0] == vIds[-1]):
                        vIds.append(dataEV[iE]["vertId"][1])
                        eIds.remove(iE)
                        break
                    elif (dataEV[iE]["vertId"][1] == vIds[-1]):
                        vIds.append(dataEV[iE]["vertId"][0])
                        eIds.remove(iE)
                        break
            vLocIds = np.array([np.argwhere(locDataV[:,2]==vId)[0,0] for vId in vIds[:-1]])
            locVerts = locDataV[vLocIds,:2]
            ax.fill(locVerts[:,0],locVerts[:,1],color=getColor(cId))
            text = str(cId)
            if (not dataCE[cId]["flat"]):
                text = text+"*"
            ax.text(cell[0],cell[1],text,color="black")
        edgeVertsCoords = {}
        def interp(scale,coords):
            x = scale*coords[1]+(1-scale)*coords[0]
            return (x[0],x[1])
        for edge in locDataE[:]:
            eId = int(edge[2])
            vLocs = [locDataV[np.argwhere(locDataV[:,2]==dataEV[eId]["vertId"][i])[0,0],:2] for i in range(len(dataEV[eId]["vertId"]))]
            edgeVertsCoords[eId] = vLocs
            ax.annotate("",xy=interp(0.52,vLocs),xytext=interp(0.48,vLocs),arrowprops={"arrowstyle": "->", "color": "black"})
        for cell in locDataC[:]:
            cId = int(cell[2])
            for iE, iO in zip(dataCE[cId]["edgeId"],dataCE[cId]["edgeOrientation"]):
                vLocs = edgeVertsCoords[iE]
                ax.annotate("",xy=interp(0.48+iO*0.24,vLocs),xytext=interp(0.48+iO*0.2,vLocs),arrowprops={"arrowstyle": "->", "color": darker(getColor(cId))})
        fig.savefig(prefix+"Orientation{}.png".format(int(mId)))
        # Plot quadrature locations
        ax.clear()
        ax.set_title('Quadrature nodes')
        for edge in locDataE[:]:
            eId = int(edge[2])
            vLocs = edgeVertsCoords[eId]
            ax.plot([vLocs[0][0],vLocs[1][0]],[vLocs[0][1],vLocs[1][1]],color="black")
        dataQuad = {}
        dataVal = {}
        with open("quadlocs.csv", "r") as fh:
            for line in fh:
                line = line.strip()
                parts = line.split(", ")
                if (mId != int(parts[0])):
                    continue
                cId = int(parts[1])
                dataQuad[cId] = np.reshape(np.array(parts[2:],dtype=np.double),(2,-1),order='F')
        with open("polyeval.csv", "r") as fh:
            for line in fh:
                line = line.strip()
                parts = line.split(", ")
                if (mId != int(parts[0])):
                    continue
                cId = int(parts[1])
                dataVal[cId] = np.reshape(np.array(parts[2:],dtype=np.double),(-1,dataQuad[cId].shape[1]),order='F')
        for cell in locDataC[:]:
            cId = int(cell[2])
            ax.plot(dataQuad[cId][0,:],dataQuad[cId][1,:],'o',color=getColor(cId))
        fig.savefig(prefix+"QuadatureLocation{}.png".format(int(mId)))
        
        polyName = []
        with open("polyname.csv", "r") as fh:
            for line in fh:
                line = line.strip()
                parts = line.split(", ")
                polyName = polyName + parts
        fig3d = plt.figure()
        fig3d.suptitle("Map Id {}".format(int(mId)))
        fig3d.set_size_inches(10,10)
        ax3d = fig3d.add_subplot(projection='3d')
        for polyBasis in range(len(polyName)):
            ax3d.set_title("Polynomial basis "+polyName[polyBasis])
            ax3d.set_xlabel('X')
            ax3d.set_ylabel('Y')
            for cell in locDataC[:]:
                cId = int(cell[2])
                ax3d.plot_trisurf(dataQuad[cId][0,:],dataQuad[cId][1,:],dataVal[cId][polyBasis,:],color=getColor(cId))
            fig3d.savefig(prefix+"quPolynomial evaluation map {} basis {}.png".format(int(mId),polyBasis))
            ax3d.clear()

if __name__ == "__main__":
    if (len(sys.argv) > 1):
        genfromcsv(sys.argv[1])

