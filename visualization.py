import ROOT as root
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cmx
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from mpl_toolkits.mplot3d import proj3d
from matplotlib.patches import FancyArrowPatch
from matplotlib import collections, colors, transforms


from detector import Detector
import physics,cutflow,util, detector




class Visualizer:
    DrawMomentum = True

    #writeDirectory = "../plots/"

    trackPtSize = 2

    def __init__(self):
        # print("Hi I'm a Visualizer")
        self.fig = plt.figure()
#        self.ax = Axes3D(self.fig)
        # self.ax = self.fig.gca(projection='3d')
        self.ax = self.fig.add_subplot(projection='3d')

        self.displayPoints = []
        self.unusedPoints = []
        self.vertices = []

        self.vert_vel = True

        self.writeDirectory = "../"


    def AddTrack(self, track):
        self.liveTracks.append(track)

    def AddPoint(self, point):
        self.displayPoints.append(point)

    def AddHit(self, point):
        self.unusedPoints.append(point)

    def AddVertex(self, point_dict):
        self.vertices.append(point_dict)

    def DetectorDisplay(self):
        det = Detector()
        xLims = det.xLims()
        zLims = det.zLims()


        self.ax.set_xlim(xLims[0], xLims[1])
        self.ax.set_ylim(det.yLims()[0], det.yLims()[1])
        self.ax.set_zlim(zLims[0], zLims[1])

        #constructing each layer for Poly3DCollection class

        layerX = [xLims[0], xLims[0], xLims[1], xLims[1]]
        layerZ = [zLims[0], zLims[1], zLims[1], zLims[0]]

        cols = []

        for layerN in range(det.numLayers()):
            layerY = [det.LayerYMid(layerN) for x in range(len(layerX))]
            verts = [list(zip(layerX, layerY, layerZ))]
            cols.append(Poly3DCollection(verts, alpha=0.10))
            cols[layerN].set_facecolor(det.DrawColor())
            self.ax.add_collection3d(cols[layerN])

        return det



    def TrackDisplay(self, ListOf_trackPtLists, Listof_colors, list_of_labels=None):
        xs, ys, zs, cs = [], [], [], []
        scatters = []
        for n, trackPtList in enumerate(ListOf_trackPtLists):
            xs += [trackPt.x for trackPt in trackPtList]
            ys += [trackPt.y for trackPt in trackPtList]
            zs += [trackPt.z for trackPt in trackPtList]
            cs += [Listof_colors[n] for trackPt in trackPtList]
            #print("cs is ",cs)
            #print("list of labels ",list_of_labels)
            scatters.append(self.ax.scatter(xs, ys, zs, s=self.trackPtSize, c=cs[n][0], label=list_of_labels[n]))
            xs, ys, zs, cs = [], [], [], []
        for pnt in self.displayPoints:
            scatters.append(self.ax.scatter(pnt[0],pnt[1],pnt[2],c=5))
        self.ax.set_xlabel("x [cm]")
        self.ax.set_ylabel("y [cm]")
        self.ax.set_zlabel("z [cm]")
        self.ax.legend(markerscale=3, loc=2)



    def TrackDisplayPoints(self, x, y, z, color=None, Label=None, opac=1):

        self.ax.plot(x, y, z, c=color, label=Label, alpha=opac)

        self.ax.set_xlabel("x [cm]")
        self.ax.set_ylabel("y [cm]")
        self.ax.set_zlabel("z [cm]")

        if Label != None:
            self.ax.legend(markerscale=3, loc=2)



    def PlotPoints(self):

        scatters = []

        for pnt in self.displayPoints:
            if len(pnt) == 3:
                scatters.append(self.ax.scatter(pnt[0],pnt[1],pnt[2],s=5,c="k",marker="*"))

            else:
                scatters.append(self.ax.scatter(pnt[0][0],pnt[0][1],pnt[0][2],s=5,c=pnt[1],marker="*"))

        for pnt in self.unusedPoints:
            if len(pnt) == 3:
                scatters.append(self.ax.scatter(pnt[0],pnt[1],pnt[2],s=10,c="k",marker="."))

            else:
                scatters.append(self.ax.scatter(pnt[0][0],pnt[0][1],pnt[0][2],s=5,c=pnt[1],marker="."))


#        print(type(self.vertices),' and ',type(dict()))
#        if type(self.vertices) == type(dict):

        #print(self.vert_vel)
        if self.vert_vel:
	        for dic in self.vertices:
	#            if len(pnt) == 3:
	             scatters.append(self.ax.scatter(dic['point'][0],dic['point'][1],dic['point'][2],s=20,c=dic['col'],marker="x"))
	             #print("with color {}".format(dic['col']))
	
	#             print(dic['vert vel'][0])
	#             print(dic['vert vel'])
	
	             if self.vert_vel:
	                   for n in range(len(dic['vert vel'][0])): # show velocity best estimates at vertex
	#                  print("velocities are {}, {}, {}".format(dic['vert vel'][0][n], dic['vert vel'][1][n], dic['vert vel'][2][n]))
	
		                   v = np.array([dic['vert vel'][0][n], dic['vert vel'][1][n], dic['vert vel'][2][n]])
		#                   print("square ",v * v)
		#                   print("initial ", v)
		                   print("Beta is ",np.sqrt(np.sum(v * v,axis=0)))
		                   v = 200 * v / np.sqrt(np.sum(v * v,axis=0))
		
		#                   print("final ",v)
		#                   print("sum ",np.sum(v * v,axis=0))
		
		
		#                                 dic['vert vel'][0][n], dic['vert vel'][1][n], dic['vert vel'][2][n],
		                   self.ax.quiver(dic['point'][0], dic['point'][1], dic['point'][2],
		                                 v[0], v[1], v[2],
		                                 color = dic['col'])# = 'k')

        else:
            for pnt in self.vertices:
                scatters.append(self.ax.scatter(pnt[0],pnt[1],pnt[2],s=20,c="r",marker="x"))

#        else:
#             v = np.array([dic['vert vel'][0], dic['vert vel'][1], dic['vert vel'][2]])
#             print("square ",v * v)
#             print("initial ", v)
#             v = 20 * v / np.sum(v * v,axis=0)
#
#             print("final ",v)
#             print("sum ",np.sum(v * v,axis=0))
#
#             self.ax.quiver(dic['point'][0], dic['point'][1], dic['point'][2],
#                                 v[0], v[1], v[2],
#                                 color = 'k')



    def Draw(self,outname='plot.pdf'):

        self.PlotPoints()

        self.DetectorDisplay()

        self.ax.view_init(elev=90,azim=-90)

        plt.savefig(self.writeDirectory + outname.split('.')[0]+'_x_.png')

        self.ax.view_init(elev=0,azim=0)

        plt.savefig(self.writeDirectory + outname.split('.')[0]+'_z_.png')


#        plt.show()





def Histogram(data, rng=None, Title=None, xaxis=None, log=False, fname='hist.png'):

    fig, ax = plt.subplots(figsize=(8,5))

    #ax.hist(data,100,range=rng)

    if log:
        ax.semilogy()


    if rng != None:
        arg_data = np.copy(data)

        data = np.array(data)

        before = np.count_nonzero(data)

        data[data < rng[0]] = 0
        data[data > rng[1]] = 0
        data[data == np.nan] = 0

        #above = len(data) - np.count_nonzero(data)
        above = before - np.count_nonzero(data) # how many got set to zero

        ax.hist(data,int(np.sqrt(len(arg_data)-above)),range=rng)

    else:
        ax.hist(data,int(np.sqrt(len(data))),range=rng)

    mean = np.mean(data)
    std = np.std(data)

    ax.set_title(Title)
    ax.set_xlabel(xaxis)

    if rng != None:
        ax.text(rng[1]*0.7,np.shape(data)[0]*5e-2,"Mean: {:.03g} \nSTD: {:0.3g} \nOverflow: {}".format(mean,std,above))

    else:
        ax.text(np.max(data)*0.7,np.shape(data)[0]*5e-2,"Mean: {:.03g} \nSTD: {:0.3g}".format(mean,std))

    plt.savefig(fname)
#    plt.show()



def root_Histogram(data, rng=None, ft_rng=None, bins=0, Title="Histogram", xaxis=None, logx=False, logy=False, fname='hist.png'):

    canv = root.TCanvas("canv","newCanvas")

    bins = int(np.sqrt(len(data))) + bins

    if rng != None:
        #hist = root.TH1F("hist",Title,bins,rng[0],rng[1])
        #hist = root.TH1F("hist",Title,bins,np.amin(data),np.max(data))
        #hist = root.TH1F("hist",Title,bins,-100,100)
        hist = root.TH1F("hist",Title+[';'+xaxis,''][xaxis==None],bins,rng[0],rng[1])

    else:
        hist = root.TH1F("hist",Title+[';'+xaxis,''][xaxis==None],bins,np.amin(data),np.max(data))

#    if len(ft_rng) != 0:
    if ft_rng != None:
        #hist.Fit('gaus','','',rng[0],rng[1])
        #fit = root.TF1("fit", "gaus", rng[0], rng[1])
        #fit = root.TF1("fit", "gaus", -3, 3)
        fit = root.TF1("fit", "gaus", ft_rng[0], ft_rng[1])

    else:
        #hist.Fit('gaus','','',np.amin(data),np.max(data))
        fit = root.TF1("fit", "gaus", np.amin(data), np.max(data))

    for elm in range(len(data)):
        hist.Fill(data[elm])

    root.gStyle.SetOptStat(111111)

    if ft_rng != None:
        hist.GetXaxis().SetRangeUser(ft_rng[0],ft_rng[1]);

        hist.Fit("fit")

        hist.GetXaxis().SetRangeUser(rng[0],rng[1]);
    
    if logx:
        canv.SetLogx()
        
    if logy:
        canv.SetLogy()

    hist.Draw()

    bin1 = hist.FindFirstBinAbove(hist.GetMaximum()/2)
    bin2 = hist.FindLastBinAbove(hist.GetMaximum()/2)
    fwhm = hist.GetBinCenter(bin2) - hist.GetBinCenter(bin1)
    hwhm = fwhm / 2

    print("fwhm is ", fwhm)
    print("hwhm is ", hwhm)

    canv.Update()
    canv.Draw()

    canv.SaveAs(fname)
    #canv.Print(fname.split('.')[0]+".pdf","PDF")



def root_2D_Histogram(data_x, data_z, xlims, zlims, Title='Plot', xbins=100, zbins=100, xlabel='x', zlabel='z', fname='plot.png'):

    if len(data_x) != len(data_z):
        return print("Length of data must be equal")

    canv = root.TCanvas("canv","newCanvas")

    hist = root.TH2F("hist",Title,xbins,xlims[0],xlims[1],zbins,zlims[0],zlims[1])

    root.gStyle.SetOptStat(111111)

#    hist.SetStats(0)
    hist.GetXaxis().SetTitle(xlabel)
    hist.GetYaxis().SetTitle(zlabel)

    for i in range(len(data_x)):
        hist.Fill(data_x[i],data_z[i])

    hist.Draw("colz")

    canv.Update()
    canv.Draw()

    canv.SaveAs(fname)
    
    
def drawdet_xz(use_cms=False, axis=None, layer_height_vis=0.2, alpha=0.1):
#     use_cms=False
#     axis=None
#     layer_height_vis=0.2
    if axis is None:
        axis=plt.gca()

    det=Detector() # Get detector geometry
    verts=[] # vertices of polygons

    # Loop 10 modules
    for ix in range(10): 
        layerX = [det.module_x_displacement[ix]-det.module_x_edge_length*0.5,
                  det.module_x_displacement[ix]-det.module_x_edge_length*0.5,
                  det.module_x_displacement[ix]+det.module_x_edge_length*0.5,
                  det.module_x_displacement[ix]+det.module_x_edge_length*0.5]
        # Loop 7 layers
        for iz in range(len(det.LayerYLims)):
            # yi = 0.5*(det.LayerYLims[i][0] + det.LayerYLims[i][1])
            # layerZ = [det.layer_z_displacement[iz]-layer_height_vis*0.5,
            #          det.layer_z_displacement[iz]+layer_height_vis*0.5,
            #          det.layer_z_displacement[iz]+layer_height_vis*0.5,
            #          det.layer_z_displacement[iz]-layer_height_vis*0.5,]
            yi = (0.5*(det.LayerYLims[iz][0] + det.LayerYLims[iz][1])-8550)*0.01 # minus the y offset and turns into meter
            layerZ = [yi-layer_height_vis*0.5,
                     yi+layer_height_vis*0.5,
                     yi+layer_height_vis*0.5,
                     yi-layer_height_vis*0.5,]            
            verts.append(np.transpose([layerX, layerZ]))

    col = collections.PolyCollection(verts, alpha=alpha)
    axis.add_collection(col)
    
drawdet_yz=drawdet_xz

def drawdet_xy(use_cms=False, axis=None, layer_height_vis=0.2, alpha=0.1):
    if axis is None:
        axis=plt.gca()

    det=Detector() # Get detector geometry
    verts=[] # vertices of polygons

    # Loop 10 modules
    for ix in range(10): 
        layerX = [det.module_x_displacement[ix]-det.module_x_edge_length*0.5,
                  det.module_x_displacement[ix]-det.module_x_edge_length*0.5,
                  det.module_x_displacement[ix]+det.module_x_edge_length*0.5,
                  det.module_x_displacement[ix]+det.module_x_edge_length*0.5]
        # Loop 10 layers
        for iy in range(10):
            layerY = [det.module_y_displacement[iy]-det.module_y_edge_length*0.5,
                      det.module_y_displacement[iy]+det.module_y_edge_length*0.5,
                      det.module_y_displacement[iy]+det.module_y_edge_length*0.5,
                      det.module_y_displacement[iy]-det.module_y_edge_length*0.5]
            verts.append(np.transpose([layerX, layerY]))

    col = collections.PolyCollection(verts, alpha=alpha)
    axis.add_collection(col)    

# Test:
# fig,ax1 = plt.subplots(1,1)
# drawdet_xz()
# plt.xlim(-50,50)
# plt.ylim(-1,12)    


cut=cutflow.sample_space("")

def plot_truth(event, fig=None, disp_det_view=True, disp_vertex=True, disp_first_hit=False):
    """
    Function to plot the truth of one event
    Tom Ren, 2023.2
    """
    # Extract truth information
    event.ExtractTruthPhysics()
    
    # Prepare the canvas
    if fig is None:
        fig,axs=plt.subplots(2,2,figsize=(12,9))
        axs=axs.flatten().tolist()
    else:
        axs=fig.axes
    
    # Plot tracks
    existing_labels=[]
    for track in event.truthTrackList:
        x = [];y = [];z = []
        for point in track.pointList:
            coord_cms = [point.location.x,point.location.y,point.location.z]
            coord_det = util.coord_cms2det(np.array(coord_cms))
            
            x.append(coord_det[0])
            y.append(coord_det[1])
            z.append(coord_det[2])
        
        # If-else to avoid duplicating lables
        # if track.LabelString() not in existing_labels:
        if 1:
            axs[0].plot(x, z, color=track.color(),marker=".",linewidth=1,markersize=4,label=track.LabelString())
            axs[1].plot(y, z, color=track.color(),marker=".",linewidth=1,markersize=4,label=track.LabelString())
            axs[2].plot(x, y, color=track.color(),marker=".",linewidth=1,markersize=4,label=track.LabelString())
            existing_labels.append(track.LabelString())
        else:
            axs[0].plot(x, z, color=track.color(),marker=".",linewidth=1,markersize=4)
            axs[1].plot(y, z, color=track.color(),marker=".",linewidth=1,markersize=4)
            axs[2].plot(x, y, color=track.color(),marker=".",linewidth=1,markersize=4)
            
    # Plot vertex
    if disp_vertex:
        try:
            used_gens_inds = np.where(np.array(event.Tree.GenParticle_G4index) != -1)[0]
            if len(used_gens_inds) != 0:
                vert_truth = [event.Tree.GenParticle_y[int(used_gens_inds[0])] / 10,
                                    event.Tree.GenParticle_x[int(used_gens_inds[0])] / 10,
                                    event.Tree.GenParticle_z[int(used_gens_inds[0])] / 10] # only first since all used

                #truth vertex location
                vertex_coord_dete = util.coord_cms2det(np.array([vert_truth[0], vert_truth[1], vert_truth[2]]))
                axs[0].scatter(vertex_coord_dete[0],vertex_coord_dete[2],s=60,marker="*", color="tab:green",alpha=0.5,zorder=100,label="Primary Vertex")
                axs[1].scatter(vertex_coord_dete[1],vertex_coord_dete[2],s=60,marker="*", color="tab:green",alpha=0.5,zorder=100)
                axs[2].scatter(vertex_coord_dete[0],vertex_coord_dete[1],s=60,marker="*", color="tab:green",alpha=0.5,zorder=100)            
        except:
            pass    
            
    if disp_first_hit:
        first_hit = util.coord_cms2det(np.array([event.Tree.Hit_x[0], event.Tree.Hit_y[0], event.Tree.Hit_z[0]]))
        axs[0].scatter(first_hit[0],first_hit[2],s=60,marker="x", color="k",alpha=0.5,zorder=100,label="First hit")
        axs[1].scatter(first_hit[1],first_hit[2],s=60,marker="x", color="k",alpha=0.5,zorder=100)
        axs[2].scatter(first_hit[0],first_hit[1],s=60,marker="x", color="k",alpha=0.5,zorder=100)

    
    if disp_det_view:
        drawdet_xz(axis=axs[0],alpha=0.2)
        drawdet_xz(axis=axs[1],alpha=0.2)
        drawdet_xy(axis=axs[2],alpha=0.2)
            
    axs[0].set_xlabel("x [m]")
    axs[0].set_ylabel("z [m]")
    axs[1].set_xlabel("y [m]")
    axs[1].set_ylabel("z [m]")
    axs[2].set_xlabel("x [m]")
    axs[2].set_ylabel("y [m]")
    # Put legend in the last grid
    handles, labels = axs[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc=(0.52,0.05))
    axs[3].axis("off")
    fig.tight_layout()
    return fig

def get_digi(event, use_cms=True):
    hits=[]
    for i in range(len(event.Tree.Digi_x)):
        hit = np.array([event.Tree.Digi_x[i],event.Tree.Digi_y[i],event.Tree.Digi_z[i]])
        hit_layer = cut.in_layer(hit[1])
        hit_uncertainty = np.array(detector.Layer().uncertainty(hit_layer))
        if not use_cms:
            hit=util.coord_cms2det(np.array(hit))
            hit_uncertainty=hit_uncertainty[[2,0,1]]
        hits.append([hit,hit_uncertainty])
    return np.array(hits)

def plot_digi(event, inds=None, fig=None, disp_det_view=False):
    """
    Function to plot the digitization of one event
    Tom Ren, 2023.2
    
    INPUT:
    inds: None or list
        indices of digitized hits to plot. Use None to plot all of them
    """
    # Extract digitized hits
    hits=get_digi(event,use_cms=False)
    
    # Prepare the canvas
    if fig is None:
        fig,axs=plt.subplots(2,2,figsize=(12,9))
        axs=axs.flatten().tolist()
    else:
        axs=fig.axes    
        
    # hits to plot:
    inds = np.arange(len(hits)) if inds is None else inds
    
    # plots hits
    x,y,z = [],[],[]
    xe,ye,ze = [],[],[]
    for i in inds:
        x.append(hits[i][0][0])
        y.append(hits[i][0][1])
        z.append(hits[i][0][2])
        xe.append(hits[i][1][0])
        ye.append(hits[i][1][1])
        ze.append(hits[i][1][2])
        
    axs[0].errorbar(x,z,xerr=xe,yerr=ze, fmt=".",capsize=2, color="red", alpha=0.3, label="digitized")
    axs[1].errorbar(y,z,xerr=ye,yerr=ze, fmt=".",capsize=2, color="red", alpha=0.3, label="digitized")
    axs[2].errorbar(x,y,xerr=xe,yerr=ye, fmt=".",capsize=2, color="red", alpha=0.3, label="digitized")
    
    if disp_det_view:
        drawdet_xz(axis=axs[0],alpha=0.2)
        drawdet_xz(axis=axs[1],alpha=0.2)
        drawdet_xy(axis=axs[2],alpha=0.2) 
    
    axs[0].set_xlabel("x [m]")
    axs[0].set_ylabel("z [m]")
    axs[1].set_xlabel("y [m]")
    axs[1].set_ylabel("z [m]")
    axs[2].set_xlabel("x [m]")
    axs[2].set_ylabel("y [m]")
    # Put legend in the last grid
    handles, labels = axs[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc=(0.52,0.05))
    axs[3].axis("off")
    fig.tight_layout()
    return fig

linestyle_tuple = [
     (0, (1, 1)),
     (0, (5, 1)),
     (0, (3, 1, 1, 1)),
     (0, (3, 1, 1, 1, 1, 1)),
     (0, (1, 4)),
     (0, (5, 5)),
     (0, (3, 5, 1, 5)),
     (0, (3, 5, 1, 5, 1, 5)),
     (0, (1, 8)),
     (5, (10, 3)),
     (0, (5, 10)),
     (0, (3, 10, 1, 10)),
     (0, (3, 10, 1, 10, 1, 10))]

def plot_recon(event, fig=None, disp_det_view=False, disp_non_vertex_tracks=True,
              disp_unused_hits=True, disp_recon_vertex=True):
    """
    Function to plot the digitization of one event
    Tom Ren, 2023.2
    
    INPUT:
    """
    # Extract reconstruction info
    event_vis = event.get_recon_kalman()
    
    # Prepare the canvas
    if fig is None:
        fig,axs=plt.subplots(2,2,figsize=(12,9))
        axs=axs.flatten().tolist()
    else:
        axs=fig.axes         
        
    key1s = ["track","track_nonvertex"]
    key2s = ["hit","hit_nonvertex"]
    if disp_non_vertex_tracks:
        loop_range = 2
    else:
        loop_range = 1
        
    # Loop both tracks used in vertex and non-vertex tracks
    for ikey in range(loop_range):
        # Plot reconstructed tracks
        name_append="" if ikey==0 else "(non vertex)"
        for i_track in range(len(event_vis[key1s[ikey]])):
            # Read the reconstructed track
            track=event_vis[key1s[ikey]][i_track]
            track=util.coord_cms2det(np.array(track))
            hits=[[],[],[]]
            hits_uncertainty=[[],[],[]]

            # Read hits of this track
            for i_hit in range(len(track[0])):
                hit=event_vis[key2s[ikey]][i_track][i_hit]
                hit_layer = cut.in_layer(hit[1])
                hit_uncertainty = np.array(detector.Layer().uncertainty(hit_layer))
                hit=util.coord_cms2det(np.array(hit))
                hit_uncertainty=hit_uncertainty[[2,0,1]]
                for i in range(3):
                    hits[i].append(hit[i])
                    hits_uncertainty[i].append(hit_uncertainty[i])
            # Plot digi hits of this track
            label = "Digitized" if i_track==0 else None
            axs[0].errorbar(hits[0],hits[2],
                                 xerr=hits_uncertainty[0],yerr=hits_uncertainty[2],
                                 color="red",capsize=2,ls='none',alpha=0.3, fmt=".", label=label)
            axs[1].errorbar(hits[1],hits[2],
                                 xerr=hits_uncertainty[1],yerr=hits_uncertainty[2],
                                 color="red",capsize=2,ls='none',alpha=0.3, fmt=".")
            axs[2].errorbar(hits[0],hits[1],
                                 xerr=hits_uncertainty[0],yerr=hits_uncertainty[1],
                                 color="red",capsize=2,ls='none',alpha=0.3, fmt=".")        
            # Plot KF reconstructed track
            axs[0].plot(track[0],track[2], color="black",linestyle=linestyle_tuple[i_track%13], linewidth=1,label=f"Recon track {i_track}{name_append}")
            axs[1].plot(track[1],track[2], color="black",linestyle=linestyle_tuple[i_track%13], linewidth=1,label=f"Recon track {i_track}{name_append}")
            axs[2].plot(track[0],track[1], color="black",linestyle=linestyle_tuple[i_track%13], linewidth=1,label=f"Recon track {i_track}{name_append}") 
            
    # Plot vertex
    if disp_recon_vertex:
        for ivertex, vertex in enumerate(event_vis["vertex"]):
            vertex_coord_dete = util.coord_cms2det(np.array(vertex))
            axs[0].scatter(vertex_coord_dete[0],vertex_coord_dete[2],s=60,marker="*", color="tab:purple",alpha=0.5,zorder=100,label=f"Recon vertex {ivertex}")
            axs[1].scatter(vertex_coord_dete[1],vertex_coord_dete[2],s=60,marker="*", color="tab:purple",alpha=0.5,zorder=100)
            axs[2].scatter(vertex_coord_dete[0],vertex_coord_dete[1],s=60,marker="*", color="tab:purple",alpha=0.5,zorder=100)  
        
    # Plot unused hits
    if disp_unused_hits:
        hits=[[],[],[]]
        hits_uncertainty=[[],[],[]]

        for i_hit in range(len(event_vis["hit_unused"])):
            hit=event_vis["hit_unused"][i_hit]
            hit_layer = cut.in_layer(hit[1])
            hit_uncertainty = np.array(detector.Layer().uncertainty(hit_layer))
            hit=util.coord_cms2det(np.array(hit))
            hit_uncertainty=hit_uncertainty[[2,0,1]]
            for i in range(3):
                hits[i].append(hit[i])
                hits_uncertainty[i].append(hit_uncertainty[i])
        axs[0].errorbar(hits[0],hits[2],xerr=hits_uncertainty[0],yerr=hits_uncertainty[2],
                             color="C0",capsize=2,ls='none',alpha=0.3, fmt=".", label="Digitized, unused" )
        axs[1].errorbar(hits[1],hits[2],xerr=hits_uncertainty[1],yerr=hits_uncertainty[2],
                             color="C0",capsize=2,ls='none',alpha=0.3, fmt=".")
        axs[2].errorbar(hits[0],hits[1],xerr=hits_uncertainty[0],yerr=hits_uncertainty[1],
                             color="C0",capsize=2,ls='none',alpha=0.3, fmt=".")      
        
    if disp_det_view:
        drawdet_xz(axis=axs[0],alpha=0.2)
        drawdet_xz(axis=axs[1],alpha=0.2)
        drawdet_xy(axis=axs[2],alpha=0.2) 
    
    axs[0].set_xlabel("x [m]")
    axs[0].set_ylabel("z [m]")
    axs[1].set_xlabel("y [m]")
    axs[1].set_ylabel("z [m]")
    axs[2].set_xlabel("x [m]")
    axs[2].set_ylabel("y [m]")
    # Put legend in the last grid
    handles, labels = axs[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc=(0.52,0.05))
    axs[3].axis("off")
    fig.tight_layout()
    return fig    
    
