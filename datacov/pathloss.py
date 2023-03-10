import numpy as np
from scipy.stats import norm
from scipy.special import lambertw
from random import choices


class MAPL:
    def __init__(self, df, n_prb):
        self.data = self.format_topo(df)
        # parameters
        self.frequency = self.data['frequency']
        self.n_prb = n_prb
        self.w = 0.180 * self.n_prb
        self.h_bs = self.data['h_bs']
        self.h_ue = 1.7
        self.T = 290
        self.kb = 1.38e-23
        self.n_d = 10*np.log10(self.kb * self.T * 1e3)
        self.p_cov = 0.95
        self.q = norm.ppf
        self.sd_sf_def = {'urban': 8, 'suburban': 7, 'rural': 6}
        self.context = ['urban', 'suburban', 'rural']
        self.set_context()
        self.set_sd_sf()
        self.sd_sf = self.data.sd_sf
        
        # transmission
        self.p_t = 23
        self.g_bs = 0
        self.eirp = self.p_t + self.g_bs
        
        # reception
        self.sinr = -6
        self.n_pow = self.n_d + 10*np.log10(self.w * 1e6)
        self.n_fig = 5
        self.g_ue_def = {
            700: 16, 
            800: 16, 
            1800: 18, 
            2100: 18, 
            2600: 19
        }
        self.set_g_ue()
        self.g_ue = self.data.g_ue
        self.l_cable = 2
        self.s = self.sinr + self.n_pow + self.n_fig - self.g_ue + self.l_cable
        
        # margins
        self.m_i = 3
        self.m_sf = self.sd_sf * self.q(self.p_cov)
        self.m_body = 1
        self.m_building_def = {'urban': 18, 'suburban': 15, 'rural': 12}
        self.set_m_building()
        self.m_building = self.data.m_building
        self.m = self.m_i + self.m_sf + self.m_body + self.m_building
        
        # mapl
        self.mapl = self.eirp - self.s - self.m
        self.set_mapl()
        #self.radius = self.rev_pl(self.mapl)
        
    def format_topo(self, topo):
        raise NotImplementedError('Not implemented.')
        
    def set_sd_sf(self):
        self.data['sd_sf'] = 0.
        for context in self.sd_sf_def.keys():
            self.data.loc[self.data['context'] == context, 'sd_sf'] = self.sd_sf_def[context]
    
    def set_context(self):
        pass
    
    def set_g_ue(self):
        self.data['g_ue'] = 0.
        for f in self.g_ue_def:
            self.data.loc[self.data.frequency == f, 'g_ue'] = self.g_ue_def[f]
            
    def set_m_building(self):
        self.data['m_building'] = 0.
        for context in self.m_building_def:
            self.data.loc[self.data['context'] == context, 'm_building'] = self.m_building_def[context]
            
    def set_mapl(self):
        self.data['mapl'] = self.mapl


class IDFMAPL(MAPL):
    def __init__(self, df, n_prb):
        super().__init__(df, n_prb)
        
    def set_context(self):
        self.data['context'] = 'suburban'
        self.data.loc[self.data['department'].isin([77, 78, 91, 95]), 'context'] = 'suburban'
        self.data.loc[self.data['TUU2017'] == 0, 'context'] = 'rural'
        self.data.loc[self.data['department'] == 75, 'context'] = 'urban'
        self.data.loc[self.data['department'].isin([92, 93, 94]), 'context'] = 'suburban'

    def format_topo(self, topo_cells):
        df = topo_cells.copy()

        aliases = {
            'hba': 'h_bs'
        }
        df['frequency'] = df['bande'].str[3:].astype(int)
        df = df.rename(columns=aliases)
        return df

class GenMAPL(MAPL):
    def __init__(self, df, n_prb, context_array):
        self.context_array = context_array
        super().__init__(df, n_prb)
        
        
    def set_context(self):
        self.data['context'] = self.context_array
        
    def format_topo(self, topo_cells):
        df = topo_cells.copy()

        aliases = {
            'height': 'h_bs',
            'band': 'frequency'
        }
        df = df.rename(columns=aliases)
        return df
        
        
class RevHata:
    def __init__(self, mapl):
        self.mapl = mapl
        self.data = self.mapl.data.copy()
        self.set_r_l()
        self.ahm = self.data.ahm
        self.r_l = self.data.r_l
    
    def set_r_l(self):
        self.data['ahm'] = 0.
        self.data['r_l'] = 0.
        for context in self.mapl.context:
            df = self.data.loc[self.data['context'] == context]
            
            if context == 'urban':
                ahm = 3.2 * np.log10(11.75 * self.mapl.h_ue)**2 - 4.97
                lgdist = (df.mapl - 69.55 - 26.16*np.log10(df.frequency) + 13.82*np.log10(df.h_bs) + ahm) / (44.9 - 6.55*np.log10(self.mapl.h_ue)) 
            if context == 'suburban':
                ahm = 3.2 * np.log10(11.75 * self.mapl.h_ue)**2 - 4.97
                lgdist = (df.mapl + 2 * np.log10(df.frequency/28)**2 + 5.4 - 69.55 - 26.16*np.log10(df.frequency) + 13.82*np.log10(df.h_bs) + ahm) / (44.9 - 6.55*np.log10(df.h_bs)) 
            if context == 'rural':
                ahm = (1.1 * np.log10(df.frequency) - 0.7) * self.mapl.h_ue - (1.56 * np.log10(df.frequency) - 0.8)
                lgdist = (df.mapl + 4.78 * np.log10(df.frequency)**2 - 18.33*np.log10(df.frequency) + 35.94 - 69.55 - 26.16*np.log10(df.frequency) + 13.82*np.log10(df.h_bs) + ahm) / (44.9 - 6.55*np.log10(df.h_bs)) 
            
            self.data.loc[self.data['context'] == context, 'ahm'] = ahm
            self.data.loc[self.data['context'] == context, 'r_l'] = 10**lgdist * 1e3
    
class RevUMa:
    def __init__(self, mapl):
        self.c = 3e8
        self.h_e = 1
        self.use_opt = True
        self.mapl = mapl
        self.data = self.mapl.data.copy()
        self.h_ut = self.mapl.h_ue
        self.h_bs = self.mapl.h_bs
        self.frequency = self.mapl.frequency
        self.set_d_bp()
        self.set_dist1()
        self.set_dist2()
        self.set_dist_los()
        self.set_dist_nlos_prime()
        self.set_dist_nlos_opt()
        self.set_dist_nlos()
        self.set_r_l()
        
    def get_d_2d(self, d_3d):
        return np.sqrt(d_3d**2 - (self.h_bs - self.h_ut)**2)
    
    def set_d_bp(self):
        self.data['d_bp'] = 4 * (self.h_bs-self.h_e).clip(0) * (self.h_ut-self.h_e) * self.frequency * 1e6 / self.c
    
    def set_dist1(self):
        # in parameters, f is asked to be in MHz
        # in formula, f is converted to GHz and d is in meters
        #context = 'urban'
        pl = self.mapl.mapl
        f = self.frequency * 1e-3
        num = pl - 28 - 20 * np.log10(f)
        dem = 22
        d_3d = 10**(num / dem)
        d_2d = self.get_d_2d(d_3d)
        self.data['dist1'] = d_2d

    def set_dist2(self):
        pl = self.mapl.mapl
        f = self.frequency * 1e-3
        num = pl - 28 - 20 * np.log10(f) + 9 * np.log10(self.data.d_bp**2 + (self.h_bs - self.h_ut)**2)
        dem = 40
        d_3d = 10**(num / dem)
        d_2d = self.get_d_2d(d_3d)
        self.data['dist2'] = d_2d
        
    def set_dist_los(self):
        d_bp = self.data['d_bp']
        d1 = self.data['dist1']
        d2 = self.data.loc[self.data['dist1'] > d_bp, 'dist2']
        
        self.data['dist_los'] = d1
        self.data.loc[self.data['dist1'] > d_bp, 'dist_los'] = d2

    def set_dist_nlos_prime(self):
        #context = 'urban'
        pl = self.mapl.mapl
        f = self.frequency * 1e-3
        num = pl - 13.54 - 20 * np.log10(f) + 0.6 * (self.h_ut-1.5)
        dem = 39.08
        d_3d = 10**(num / dem)
        d_2d = self.get_d_2d(d_3d)
        self.data['dist_nlos_prime'] = d_2d

    def set_dist_nlos_opt(self):
        pl = self.mapl.mapl
        f = self.frequency * 1e-3
        num = pl - 32.4 - 20 * np.log10(f)
        dem = 30
        d_3d = 10**(num / dem)
        d_2d = self.get_d_2d(d_3d)
        self.data['dist_nlos_opt'] = d_2d

    def set_dist_nlos(self):
        if self.use_opt is True:
            self.data['dist_nlos'] = self.data['dist_nlos_opt']
        else:
            d1 = self.data['dist_nlos_prime']
            d2 = self.data.loc[self.data['dist_nlos_prime'] > 5000, 'dist_nlos_opt']

            self.data['dist_nlos'] = d1
            self.data.loc[self.data['dist_nlos_prime'] > 5000, 'dist_nlos'] = d2
        
    def set_r_l(self):
        if self.use_opt is True:
            self.data['r_l'] = self.data['dist_los']
            dist_prime = self.data.loc[self.data['context'].isin(['urban']), 'dist_nlos_prime']
            dist_opt = self.data.loc[self.data['context'].isin(['suburban']), 'dist_nlos_opt']
            
            self.data.loc[self.data['context'].isin(['urban']), 'r_l'] = dist_prime
            self.data.loc[self.data['context'].isin(['suburban']), 'r_l'] = dist_opt
            
        else:
            self.data['r_l'] = self.data['dist_los']
            d = self.data.loc[self.data['context'].isin(['urban', 'suburban']), 'dist_nlos']
            self.data.loc[self.data['context'].isin(['urban', 'suburban']), 'r_l'] = d

    # what we want is actually line of sight scenario
    
class RevRMa:
    def __init__(self, mapl):
        self.c = 3e8
        self.h_e = 1
        self.h = 5 # avg building height
        self.street_width = 10
        self.mapl = mapl
        self.data = self.mapl.data.copy()
        self.h_ut = self.mapl.h_ue
        self.h_bs = self.mapl.h_bs
        self.frequency = self.mapl.frequency
        
        self.set_d_bp()
        self.set_dist1()
        self.set_dist2()
        self.set_dist_los()
        self.set_dist_nlos_prime()
        self.set_r_l()
        
    def get_d_2d(self, d_3d):
        return np.sqrt(d_3d**2 - (self.h_bs - self.h_ut)**2)
    
    def set_d_bp(self):
        self.data['d_bp'] = 2 * np.pi * self.h_bs * self.h_ut * self.frequency * 1e6 / self.c
    
    def set_dist1(self):
        # in parameters, f is asked to be in MHz
        # in formula, f is converted to GHz and d is in meters
        #context = 'urban'
        pl = self.mapl.mapl
        f = self.frequency * 1e-3
        
        # with RMa, PL function is in the form Alog(x) + Bx + C
        # we transform it to be in a form (PL-K)/B = A/Blog(x) + x
        K = 20*np.log10(40*np.pi*f/3) - 0.44*self.h**1.72
        A = 20 + 0.03*self.h**1.72
        B = 0.002*np.log10(self.h)
        
        X = (B / A * np.log(10) * 10**((pl-K)/A)).to_numpy()
        w_X = lambertw(X).real
        d_3d = (A / (B * np.log(10)) * w_X)
        
        d_2d = self.get_d_2d(d_3d)
        self.data['dist1'] = d_2d

    def set_dist2(self):
        pl = self.mapl.mapl
        f = self.frequency * 1e-3
        d_bp = self.data.d_bp
        
        K = 20*np.log10(40*np.pi*f/3) - 0.044*self.h**1.72 - 40 * np.log10(d_bp)
        A = 60+0.03*self.h**1.72
        B = 0.002*np.log10(self.h)
        
        
        X = (B / A * np.log(10) * 10**((pl-K)/A)).to_numpy()
        w_X = lambertw(X).real
        d_3d = (A / (B * np.log(10)) * w_X)
        
        d_2d = self.get_d_2d(d_3d)
        self.data['dist2'] = d_2d
        
    def set_dist_los(self):
        d_bp = self.data['d_bp']
        d1 = self.data['dist1']
        d2 = self.data.loc[self.data['dist1'] > d_bp, 'dist2']
        
        self.data['dist_los'] = d1
        self.data.loc[self.data['dist1'] > d_bp, 'dist_los'] = d2

    def set_dist_nlos_prime(self):
        #context = 'urban'
        pl = self.mapl.mapl
        f = self.frequency * 1e-3
        
        num = pl - (161.04 - 7.1*np.log10(self.street_width) + 7.5*np.log10(self.h) - (24.37 - 3.7 * (self.h/self.h_bs)**2) * np.log10(self.h_bs) - 3 * (43.42 - 3.1*np.log10(self.h_bs)) + 20*np.log10(f) - (3.2 * (np.log10(11.75*self.h_ut))**2 - 4.97))
        dem = 43.42-3.1*np.log10(self.h_bs)
        d_3d = 10**(num / dem)
        d_2d = self.get_d_2d(d_3d)
        self.data['dist_nlos_prime'] = d_2d
        
    def set_r_l(self):
        self.data['r_l'] = self.data['dist_los']
        dist_prime = self.data.loc[self.data['context'].isin(['urban', 'suburban']), 'dist_nlos_prime']
        self.data.loc[self.data['context'].isin(['urban', 'suburban']), 'r_l'] = dist_prime
        

class RevMa:
    def __init__(self, mapl):
        self.uma = RevUMa(mapl)
        self.rma = RevRMa(mapl)
        self.data = mapl.data.copy()
        self.set_r_l()
        
    def set_r_l(self):
        self.data['r_l'] = self.rma.data['dist_los']
        d = self.uma.data.loc[self.data['context'].isin(['urban', 'suburban']), 'dist_nlos']
        self.data.loc[self.data['context'].isin(['urban', 'suburban']), 'r_l'] = d