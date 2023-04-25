from tqdm import tqdm
import iminuit
from iminuit import Minuit, cost
import numpy as np
import scipy as sp


import cutflow,detector,event,util
cut=cutflow.sample_space("")
det=detector.Detector()

class Hit:
    def __init__(self, x, y, z, t):
        self.x=x
        self.y=y
        self.z=z
        self.t=t
        self.t_uncertainty=1
    def get_uncertainty(self):
        # Get the layer-dependent uncertainty of each hit
        self.hit_layer=cut.in_layer(self.y)
        hit_uncertainty = np.array(detector.Layer().uncertainty(self.hit_layer))
        self.x_uncertainty=hit_uncertainty[0]*100 # m->cm
        self.z_uncertainty=hit_uncertainty[2]*100 # m->cm
        self.y_uncertainty=2/np.sqrt(12)

def get_digi_hits(ev):
    ev.Tree.GetEntry(ev.EventNumber)
    hits=[]
    for ii in range(len(ev.Tree.Digi_y)):
        hit=Hit(ev.Tree.Digi_x[ii], ev.Tree.Digi_y[ii], ev.Tree.Digi_z[ii], ev.Tree.Digi_time[ii])
        hit.get_uncertainty()
        hits.append(hit)
    return hits

def get_event_truth(ev):
    ev.Tree.GetEntry(ev.EventNumber)
    dx=ev.Tree.Hit_x[1]-ev.Tree.Hit_x[0]
    dy=ev.Tree.Hit_y[1]-ev.Tree.Hit_y[0]
    dz=ev.Tree.Hit_z[1]-ev.Tree.Hit_z[0]
    dt=ev.Tree.Hit_time[1]-ev.Tree.Hit_time[0]
    truth=[ev.Tree.Hit_x[0], ev.Tree.Hit_y[0], ev.Tree.Hit_z[0], ev.Tree.Hit_time[0],dx/dt, dy/dt, dz/dt]
    return truth




# -------------------------------------
# LS fit
# ------------------------------------
class chi2_track:
    def __init__(self, hits):
        self.hits=hits
        self.func_code = iminuit.util.make_func_code(['x0', 'y0', 'z0', 't0', 'vx', 'vy', 'vz'])
    def __call__(self, x0, y0, z0, t0, vx, vy, vz):
        error=0
        for hit in self.hits:
            model_t = (hit.y - y0)/vy
            model_x = x0 + model_t*vx
            model_z = z0 + model_t*vz
            error+= np.sum(np.power([(model_t- (hit.t-t0))/hit.t_uncertainty, 
                                     (model_x-hit.x)/hit.x_uncertainty, 
                                     (model_z-hit.z)/hit.z_uncertainty],2))
        return error        

def guess_track(hits):
    # Guess initial value
    x0_init = hits[0].x
    y0_init = hits[0].y
    z0_init = hits[0].z
    t0_init = hits[0].t
    dt=hits[-1].t-hits[0].t
    vx_init = (hits[-1].x-hits[0].x)/dt
    vy_init = (hits[-1].y-hits[0].y)/dt
    vz_init = (hits[-1].z-hits[0].z)/dt
    v_mod = np.sqrt(vx_init**2+vy_init**2+vz_init**2)
    if v_mod>sp.constants.c*1e-7:
        vx_init = vx_init*0.99*sp.constants.c*1e-7/v_mod
        vy_init = vy_init*0.99*sp.constants.c*1e-7/v_mod
        vz_init = vz_init*0.99*sp.constants.c*1e-7/v_mod
    return  (x0_init,y0_init, z0_init,t0_init,vx_init,vy_init,vz_init)
    
def fit_track(hits, guess):
    x0_init,y0_init, z0_init,t0_init,vx_init,vy_init,vz_init = guess
    det=detector.Detector()

    m = Minuit(chi2_track(hits),x0=x0_init, y0=y0_init, z0=z0_init, t0=t0_init, vx=vx_init, vy=vy_init, vz=vz_init)
    m.fixed["y0"]=True
    m.limits["x0"]=(det.BoxLimits[0][0],det.BoxLimits[0][1])
    m.limits["z0"]=(det.BoxLimits[2][0],det.BoxLimits[2][1])
    m.limits["t0"]=(-100,1e5)
    m.limits["vx"]=(-sp.constants.c*1e-7, sp.constants.c*1e-7) # Other
    m.limits["vy"]=(-sp.constants.c*1e-7*2,0) if vy_init<0 else (0,sp.constants.c*1e-7*2) # Constrain the direction in Z(up) in real world
    m.limits["vz"]=(-sp.constants.c*1e-7, sp.constants.c*1e-7) # Beam direction; From MKS unit to cm/ns = 1e2/1e9=1e-7
    m.errors["x0"]=0.1
    m.errors["y0"]=0.1
    m.errors["z0"]=0.1
    m.errors["t0"]=0.3
    m.errors["vx"] = 0.01
    m.errors["vy"] = 0.1
    m.errors["vz"] = 0.01

    m.migrad()  # run optimiser
    m.hesse()   # run covariance estimator
    
    return m

def fit_track_scipy(hits, guess):
    # Using scipy?
    bounds=[]
    bounds.append((det.BoxLimits[0][0],det.BoxLimits[0][1]))
    bounds.append((det.BoxLimits[2][0],det.BoxLimits[2][1]))
    bounds.append((0,1e10))
    bounds.append((-sp.constants.c*1e-7, sp.constants.c*1e-7))
    lim_vy=(-sp.constants.c*1e-7,0) if guess[-2]<0 else (0,sp.constants.c*1e-7) # Constrain the direction in Z(up) in real world
    bounds.append(lim_vy)
    bounds.append((-sp.constants.c*1e-7, sp.constants.c*1e-7)) # Beam direction; From MKS unit to cm/ns = 1e2/1e9=1e-7    
    
    def ls(par): 
        return chi2_track(hits)(par[0], guess[1], *par[1:])
    guess_no_y=guess[:1]+guess[2:]
    res = sp.optimize.minimize(ls, x0=guess_no_y,method="powell",bounds=bounds)  
    
    return res


def do_ls(filename, nfit=None):
    # tfile = root.TFile.Open(filename)
    # Tree = tfile.Get("integral_tree")
    
    ev = event.Event(filename, 0, tree_name="integral_tree")
    Tree=ev.Tree
    cut=cutflow.sample_space("")
    nevents = int(Tree.GetEntries())    

    truth=[]
    truth_nlayer=[]
    
    recon=[]
    recon_unc=[]
    recon_ls=[]
    recon_ls_unc=[]
    
    ndigi_total=[]
    ndigi_track=[]
    ndigi_inds=[]
    ntracks=[]
    mask_recon_success=[]
    
    if nfit is None:
        nfit=nevents
    for Entry in tqdm(range(nfit)):
    # for Entry in range(40):
        #Tree.GetEntry(Entry)
        ev.EventNumber=Entry
        ev.Tree.GetEntry(Entry)
        hits = get_digi_hits(ev)
                
        # Get truth (speed need to be calculated by hand)
        try:
            # Truth position and speed
            dt=Tree.Hit_time[1]-Tree.Hit_time[0]
            vx=(Tree.Hit_x[1]-Tree.Hit_x[0])/dt
            vy=(Tree.Hit_y[1]-Tree.Hit_y[0])/dt
            vz=(Tree.Hit_z[1]-Tree.Hit_z[0])/dt
            truth.append([Tree.Hit_z[0], Tree.Hit_x[0], Tree.Hit_y[0], Tree.Hit_time[0],vz,vx,vy])  
            
            # Truth number of layer the first partical goes through
            pdgids = np.array([Tree.Hit_particlePdgId[i] for i in range(len(Tree.Hit_particlePdgId))])
            ind_lasthit = int(np.argmax(pdgids!=pdgids[0]))-1
            y_lasthit = Tree.Hit_y[ind_lasthit]
            y_layer = cut.in_layer(y_lasthit)
            truth_nlayer.append(y_layer)      
            
        except:
            truth.append([-9999, -9999, -9999, -9999, -9999, -9999, -9999])
            truth_nlayer.append(-9999)
           
            
        
        # If there is reconstruction:
        if len(Tree.Track_k_m_z0)==0:
            recon.append([-9990, -9990, -9990, -9990, -9990, -9990, -9990])
            recon_unc.append([-9990, -9990, -9990, -9990, -9990, -9990, -9990])
            recon_ls.append([-9990, -9990, -9990, -9990, -9990, -9990, -9990])
            recon_ls_unc.append([-9990, -9990, -9990, -9990, -9990, -9990, -9990])            
            ndigi_total.append(-9999)
            ndigi_track.append(-9999)
            ndigi_inds.append([-9999])
            ntracks.append(-9999)
            mask_recon_success.append(False)
        else:
            # Check which one is closest to truth
            track_digi_hit_inds = util.unzip(Tree.Track_k_m_hitIndices)
            track_digi_hit_len = np.array([len(i) for i in track_digi_hit_inds])
            track_chi2s = []
            
            if len(track_digi_hit_inds)>1:
                for track_ind in range(len(track_digi_hit_inds)):
                    recon_i = [Tree.Track_k_m_z0[track_ind], Tree.Track_k_m_x0[track_ind], Tree.Track_k_m_y0[track_ind], Tree.Track_k_m_t0[track_ind],Tree.Track_k_m_velZ[track_ind], Tree.Track_k_m_velX[track_ind], Tree.Track_k_m_velY[track_ind]]
                    recon_i_unc = [Tree.Track_k_m_ErrorZ0[track_ind], Tree.Track_k_m_ErrorX0[track_ind], Tree.Track_k_m_ErrorY0[track_ind], Tree.Track_k_m_ErrorT0[track_ind],Tree.Track_k_m_ErrorVz[track_ind], Tree.Track_k_m_ErrorVx[track_ind], Tree.Track_k_m_ErrorVy[track_ind]]
                    chi2 = util.chi2_calc(recon_i,truth[-1],recon_i_unc)
                    track_chi2s.append(chi2)
                #     print("s")
                #     print(track_chi2s)
                track_ind = int(np.argmin(track_chi2s))
            else:
                track_ind=0
            recon_i = [Tree.Track_k_m_z0[track_ind], Tree.Track_k_m_x0[track_ind], Tree.Track_k_m_y0[track_ind], Tree.Track_k_m_t0[track_ind],Tree.Track_k_m_velZ[track_ind], Tree.Track_k_m_velX[track_ind], Tree.Track_k_m_velY[track_ind]]
            recon_i_unc = [Tree.Track_k_m_ErrorZ0[track_ind], Tree.Track_k_m_ErrorX0[track_ind], Tree.Track_k_m_ErrorY0[track_ind], Tree.Track_k_m_ErrorT0[track_ind],Tree.Track_k_m_ErrorVz[track_ind], Tree.Track_k_m_ErrorVx[track_ind], Tree.Track_k_m_ErrorVy[track_ind]]
            
            recon.append(recon_i)
            recon_unc.append(recon_i_unc)
            ndigi_total.append(len(Tree.Digi_x))
            ndigi_track.append(track_digi_hit_len[track_ind])
            ndigi_inds.append(track_digi_hit_inds[track_ind])
            ntracks.append(len(track_digi_hit_inds))
            mask_recon_success.append(True)                
            
            
            # Do LS fit
            hits_fit=np.array(hits)[track_digi_hit_inds[track_ind]]
            guess=guess_track(hits_fit)
            fit1=fit_track(hits_fit,guess)
            par_fit=np.array(list(fit1.values))
            par_fit_error=np.array(list(fit1.errors))
            # Save results
            recon_ls.append(par_fit[[2,0,1,3,6,4,5]])
            recon_ls_unc.append(par_fit_error[[2,0,1,3,6,4,5]])
            
                
    results={
        "truth":np.array(truth),
        "truth_nlayer":np.array(truth_nlayer),
        "recon":np.array(recon),
        "recon_unc":np.array(recon_unc),
        "recon_ls":np.array(recon_ls),
        "recon_ls_unc":np.array(recon_ls_unc),        
        "ndigi_total":np.array(ndigi_total),
        "ndigi_track":np.array(ndigi_track),
        "ndigi_inds":np.array(ndigi_inds),
        "mask_recon_success":np.array(mask_recon_success)
    }
    
    return results



def get_km(filename, results_fit, tree_name="integral_tree"):
    results_fit["par_km"]=[]
    results_fit["par_km_error"]=[]
    results_fit["par_km_truth"]=[]
    results_fit["par_km_ndigi"]=[]
    results_fit["par_km_ndigitrack"]=[]
    results_fit["par_km_chi2"]=[]
    
    ev = event.Event(filename, 0, tree_name=tree_name)
    nevents = int(ev.Tree.GetEntries())
    cut=cutflow.sample_space("")


    for i_event in tqdm(range(nevents)):
        ev.EventNumber=i_event
        ev.Tree.GetEntry(i_event)

        hits = get_digi_hits(ev)
        nhits=len(ev.Tree.Hit_x)
        
        # Use Try to only process events with kalman reconstruction
        try:
            digi_hit_inds = util.unzip(ev.Tree.Track_k_m_hitIndices)
            digi_hit_len = np.array([len(i) for i in digi_hit_inds])
            track_ind = int(np.argmax(digi_hit_len==7))
            track_hits_inds=digi_hit_inds[track_ind]   
            hits_fit=np.array(hits)[track_hits_inds]
        
            par_km =[ev.Tree.Track_k_m_x0[track_ind], ev.Tree.Track_k_m_y0[track_ind], ev.Tree.Track_k_m_z0[track_ind], ev.Tree.Track_k_m_t0[track_ind], ev.Tree.Track_k_m_velX[track_ind], ev.Tree.Track_k_m_velY[track_ind], ev.Tree.Track_k_m_velZ[track_ind]]
            par_km_error =[ev.Tree.Track_k_m_ErrorX0[track_ind], ev.Tree.Track_k_m_ErrorY0[track_ind], ev.Tree.Track_k_m_ErrorZ0[track_ind], ev.Tree.Track_k_m_ErrorT0[track_ind], ev.Tree.Track_k_m_ErrorVx[track_ind], ev.Tree.Track_k_m_ErrorVy[track_ind], ev.Tree.Track_k_m_ErrorVz[track_ind]]
            xyz0_km = util.coord_cms2det(np.array(par_km[:3]))
            xyzV_km = np.array([par_km_error[2+4],par_km_error[0+4],-par_km_error[1+4]])
            
        except:
            continue    
        results_fit["par_km"].append(par_km)
        results_fit["par_km_error"].append(par_km_error)
        results_fit["par_km_ndigi"].append(len(ev.Tree.Digi_x))
        results_fit["par_km_ndigitrack"].append(digi_hit_len[track_ind])
        results_fit["par_km_chi2"].append(ev.Tree.Track_k_m_smooth_chi_sum[track_ind])
        
        Tree=ev.Tree
        dt=Tree.Hit_time[1]-Tree.Hit_time[0]
        vx=(Tree.Hit_x[1]-Tree.Hit_x[0])/dt
        vy=(Tree.Hit_y[1]-Tree.Hit_y[0])/dt
        vz=(Tree.Hit_z[1]-Tree.Hit_z[0])/dt
        par_truth = [Tree.Hit_x[0], Tree.Hit_y[0], Tree.Hit_z[0], Tree.Hit_time[0],vx,vy,vz]
        results_fit["par_km_truth"].append(par_truth)        
        
    # results_fit["par_km"]=np.array(results_fit["par_km"])
    # results_fit["par_km_error"]=np.array(results_fit["par_km_error"])
    for key in results_fit:
        results_fit[key]=np.array(results_fit[key])