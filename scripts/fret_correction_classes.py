import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd


"""
WARNING:
- All spectra used for this analysis should be previously scaled and interpolated so that they are all at the same resolution (same step-width) and within the same range (min and max)

"""
class fluorophore:
    """
    Definiton a fluorophore containing all the optical information needed for corrections
    name = simply the name of the fluorophore
    qy = measured or reported quantum yield for the fluorophore
    espilon = molar excitinction coefficient for the fluorophore
    exc_spec = excitation spectrum of the fluorophore
    emi_spec = emission spectrum of the fluorophore
    """
    def __init__(self, name, qy, epsilon, exc_spec, emi_spec, wref):
        self.name = name # name of fluorophore without SNAP/HALO tag
        self.qy = float(qy) # quantum yield of fluorophore
        self.epsilon = float(epsilon) # molar exctinction coefficient of fluorophore
        self.exc_spec = exc_spec # excitation interpolated spectrum
        self.emi_spec = emi_spec # emission interpolated spectrum
        self.wref = wref.to_numpy() # reference wavelenghts for which spectra are interpolated to
        
        # normalize emission spectra so area under the curve = QY
        self.emi_spec_norm = (self.emi_spec / np.trapz(self.emi_spec)) * self.qy

class couple:
    """
    Definition of a donor-acceptor couple for FRET corrections
    Create object with coupleX = couple(params..)
    Check Excitation efficiency with coupleX.plot_excitation()
    Check Emission efficiency with coupleX.plot_emission()
    coupleX.calc_correction_factors(plot=True/False) will print out calculated correction parameters with optional plotting of how those parameters are obtained
    """
    def __init__(self, couple_name, dd_raw, da_raw, aa_raw, d, a, d_dm, a_dm, d_emi_filter, a_emi_filter, d_led, a_led, d_exc_filter, a_exc_filter, d_pow, a_pow, camera, dd_exposure, da_exposure, aa_exposure, wref, d_dm_ref=None, a_dm_ref=None):
        self.name = couple_name # string to keep as name reference. Provide format to be displayed in plots
        self.dd_raw = dd_raw # pandas dataframe containing traces of donor excitation, donor emission
        self.da_raw = da_raw # pandas dataframe containing traces of donor excitation, acceptor emission
        self.aa_raw = aa_raw # pandas dataframe containing traces of acceptor excitation, acceptor emission
        self.d = d # donor fluorophore object
        self.a = a # acceptor fluorophore object
        self.d_dm = d_dm # dichlabelc mirror used for donor
        self.a_dm = a_dm # dichlabelc mirror used for acceptor
        self.d_emi_filter = d_emi_filter # emission filter used for donor
        self.a_emi_filter = a_emi_filter # emission filter used for acceptor
        self.d_led = d_led # ligth source spectrum used for donor excitation
        self.a_led = a_led # ligth source spectrum used for acceptor excitation
        self.d_exc_filter = d_exc_filter # excitation filter used for donor
        self.a_exc_filter = a_exc_filter # excitation filter used for acceptor
        self.d_pow = d_pow # ligth source power used for donor excitation with corresponding DM and filters
        self.a_pow = a_pow # ligth source power used for acceptor excitation with corresponding DM and filters
        self.camera = camera # camera detection efficiency spectrum
        self.dd_exposure = dd_exposure
        self.da_exposure = da_exposure
        self.aa_exposure = aa_exposure
        self.wref = wref.to_numpy() # reference wavelenghts for which spectra are interpolated to
        
        # automatically calculate the reflexion of dichlabelc mirrors
        self.d_dm_ref = d_dm_ref if d_dm_ref is not None else np.ones(self.wref.shape[0], dtype='float') - self.d_dm
        self.a_dm_ref = a_dm_ref if a_dm_ref is not None else np.ones(self.wref.shape[0], dtype='float') - self.a_dm

    def plot_excitation(self):
        fig, axes = plt.subplots(3,2, figsize=(14,8))
        
        axes[0,0].set_title('DONOR channel')
        axes[0,0].plot(self.wref, self.d_dm, label=self.d_dm.name)
        axes[0,0].plot(self.wref, self.d_exc_filter, label=self.d_exc_filter.name)
        axes[0,0].plot(self.wref, self.d_led, label=self.d_led.name)
        axes[0,0].plot(self.wref, self.d.exc_spec, label='excitation {}'.format(self.d.name))
        
        d_max_pos = self.d.exc_spec.argmax()
        a_max_pos = self.a.exc_spec.argmax()
        d_max_wl = self.wref[d_max_pos]
        a_max_wl = self.wref[a_max_pos]
        
        axes[1,0].plot(self.wref, self.d_dm, label=self.d_dm.name)
        axes[1,0].plot(self.wref, self.d_exc_filter, label=self.d_exc_filter.name)
        axes[1,0].plot(self.wref, self.d_led, label=self.d_led.name)
        axes[1,0].plot(self.wref, self.a.exc_spec, label='excitation {}'.format(self.a.name))

        DDd_effective_excitation = self.d_dm_ref * self.d_exc_filter * self.d_led * self.d.exc_spec
        DDa_effective_excitation = self.d_dm_ref * self.d_exc_filter * self.d_led * self.a.exc_spec
        axes[2,0].plot(self.wref, DDd_effective_excitation, label='DDd effective excitation')
        axes[2,0].plot(self.wref, DDa_effective_excitation, label='DDa effective excitation')
        
        axes[0,1].set_title('ACCEPTOR channel')
        axes[0,1].plot(self.wref, self.a_dm, label=self.a_dm.name)
        axes[0,1].plot(self.wref, self.a_exc_filter, label=self.a_exc_filter.name)
        axes[0,1].plot(self.wref, self.a_led, label=self.a_led.name)
        axes[0,1].plot(self.wref, self.a.exc_spec, label='excitation {}'.format(self.a.name))        
        
        axes[1,1].plot(self.wref, self.a_dm, label=self.a_dm.name)
        axes[1,1].plot(self.wref, self.a_exc_filter, label=self.a_exc_filter.name)
        axes[1,1].plot(self.wref, self.a_led, label=self.a_led.name)
        axes[1,1].plot(self.wref, self.d.exc_spec, label='excitation {}'.format(self.d.name)) 

        AAa_effective_excitation = self.a_dm_ref * self.a_exc_filter * self.a_led * self.a.exc_spec
        AAd_effective_excitation = self.a_dm_ref * self.a_exc_filter * self.a_led * self.d.exc_spec
        axes[2,1].plot(self.wref, AAa_effective_excitation, label='AAa effective excitation')
        axes[2,1].plot(self.wref, AAd_effective_excitation, label='AAd effective excitation')
        
        for ax in axes.flat:
            ax.set_xlim(d_max_wl-100, a_max_wl+100)
            ax.legend()
            
    def plot_emission(self):

        fig, axes = plt.subplots(3,2, figsize=(14,8))

        axes[0,0].set_title('DONOR channel')
        axes[0,0].plot(self.wref, self.d_dm, label=self.d_dm.name)
        axes[0,0].plot(self.wref, self.d_emi_filter, label=self.d_emi_filter.name)
        axes[0,0].plot(self.wref, self.d.emi_spec, label='emission {}'.format(self.d.name))
        
        axes[1,0].plot(self.wref, self.d_dm, label=self.d_dm.name)
        axes[1,0].plot(self.wref, self.d_emi_filter, label=self.d_emi_filter.name)
        axes[1,0].plot(self.wref, self.a.emi_spec, label='emission {}'.format(self.a.name))

        d_max_pos = self.d.emi_spec.argmax()
        a_max_pos = self.a.emi_spec.argmax()
        self.d_max_wl = self.wref[d_max_pos]
        self.a_max_wl = self.wref[a_max_pos]
        
        DDd_effective_emission = self.d_dm * self.d_emi_filter * self.d.emi_spec
        DDa_effective_emission = self.d_dm * self.d_emi_filter * self.a.emi_spec
        axes[2,0].plot(self.wref, DDd_effective_emission, label='DDd effective emission')
        axes[2,0].plot(self.wref, DDa_effective_emission, label='DDa effective emission')
        
        axes[0,1].set_title('ACCEPTOR channel')
        axes[0,1].plot(self.wref, self.a_dm, label=self.a_dm.name)
        axes[0,1].plot(self.wref, self.a_emi_filter, label=self.a_emi_filter.name)
        axes[0,1].plot(self.wref, self.a.emi_spec, label='emission {}'.format(self.a.name))
        
        axes[1,1].plot(self.wref, self.a_dm, label=self.a_dm.name)
        axes[1,1].plot(self.wref, self.a_emi_filter, label=self.a_emi_filter.name)
        axes[1,1].plot(self.wref, self.d.emi_spec, label='emission {}'.format(self.d.name))

        AAa_effective_emission = self.a_dm * self.a_emi_filter * self.a.emi_spec
        AAd_effective_emission = self.a_dm * self.a_emi_filter * self.d.emi_spec
        axes[2,1].plot(self.wref, AAa_effective_emission, label='AAa effective emission')
        axes[2,1].plot(self.wref, AAd_effective_emission, label='AAd effective emission')
        
        for ax in axes.flat:
            ax.set_xlim(self.d_max_wl-100, self.a_max_wl+100)
            ax.legend()

    def calc_correction_factors_full(self, quiet=False, plot=False):

        ### Normalize to measured light source power
        # We take  into account all filters through which the power was  measured
        # With self.wref[DD_filters.argmax()] factor we are correcting for hv the photons amount for wavelength
        DD_filters = self.d_dm_ref * self.d_led * self.d_exc_filter
        DD_filters_norm = (self.d_pow * DD_filters * self.wref[DD_filters.argmax()]) / np.trapz(DD_filters)

        AA_filters = self.a_dm_ref * self.a_led * self.a_exc_filter
        AA_filters_norm = (self.a_pow * AA_filters * self.wref[DD_filters.argmax()]) / np.trapz(AA_filters)
        
        ## Sdd components      
        self.Ex_DDd = self.d.epsilon * DD_filters_norm * self.d.exc_spec # how much donor is excited for DONOR channel?
        self.Ex_DDd_graph = self.d_dm_ref * self.d_led * self.d_exc_filter * self.d.exc_spec
        self.Em_DDd = self.d_dm * self.d_emi_filter * self.d.emi_spec_norm * self.camera # how much donor is in DONOR channel?
        self.Em_DDd_graph = self.d_dm * self.d_emi_filter * self.d.emi_spec * self.camera
        self.Ex_DDa = self.a.epsilon * DD_filters_norm * self.a.exc_spec # how much acceptor is excited for DONOR channel?
        self.Ex_DDa_graph = self.d_dm_ref * self.d_led * self.d_exc_filter * self.a.exc_spec
        self.Em_DDa = self.d_dm * self.d_emi_filter * self.a.emi_spec_norm * self.camera # how much acceptor is Donor channel?
        self.Em_DDa_graph = self.d_dm * self.d_emi_filter * self.a.emi_spec * self.camera

        ## Sda components
        self.Em_DAd = self.d_dm * self.a_emi_filter * self.d.emi_spec_norm * self.camera
        self.Em_DAd_graph = self.d_dm * self.a_emi_filter * self.d.emi_spec * self.camera # how much donor is in FRET channel?
        self.Em_DAa = self.d_dm * self.a_emi_filter * self.a.emi_spec_norm * self.camera # how much acceptor is in FRET channel?
        self.Em_DAa_graph = self.d_dm * self.a_emi_filter * self.a.emi_spec * self.camera

        ## Saa components

        self.Ex_AAd = self.d.epsilon * AA_filters_norm * self.d.exc_spec # how much donor is excited for ACCEPTOR channel?
        self.Ex_AAd_graph = self.a_dm_ref * self.a_led * self.a_exc_filter * self.d.exc_spec
        self.Em_AAd = self.a_dm * self.a_emi_filter * self.d.emi_spec_norm * self.camera # how much donor is in ACCEPTOR channel?
        self.Em_AAd_graph = self.a_dm * self.a_emi_filter * self.d.emi_spec * self.camera
        self.Em_AAa = self.a_dm * self.a_emi_filter * self.a.emi_spec_norm * self.camera # how much acceptor is in ACCEPTOR channel?
        self.Em_AAa_graph = self.a_dm * self.a_emi_filter * self.a.emi_spec * self.camera
        self.Ex_AAa = self.a.epsilon * AA_filters_norm * self.a.exc_spec # how much acceptor is excited for ACCEPTOR channel?
        self.Ex_AAa_graph = self.a_dm_ref * self.a_led * self.a_exc_filter * self.a.exc_spec

        ## Conditional to prevent x parameters to be zero. Replace with very low number instead
        self.x0 = np.trapz(self.Ex_DDd) * np.trapz(self.Em_DDd) if np.trapz(self.Ex_DDd) * np.trapz(self.Em_DDd) else 1e-10
        self.x1 = (np.trapz(self.Em_DDa) - np.trapz(self.Em_DDd)) * np.trapz(self.Ex_DDd) if (np.trapz(self.Em_DDa) - np.trapz(self.Em_DDd)) * np.trapz(self.Ex_DDd) else 1e-10
        self.x2 = np.trapz(self.Ex_DDa) * np.trapz(self.Em_DDa) if np.trapz(self.Ex_DDa) * np.trapz(self.Em_DDa) else 1e-10
        self.x3 = np.trapz(self.Ex_DDd) * np.trapz(self.Em_DAd) if np.trapz(self.Ex_DDd) * np.trapz(self.Em_DAd) else 1e-10
        self.x4 = (np.trapz(self.Em_DAa) - np.trapz(self.Em_DAd)) * np.trapz(self.Ex_DDd) if (np.trapz(self.Em_DAa) - np.trapz(self.Em_DAd)) * np.trapz(self.Ex_DDd) else 1e-10
        self.x5 = np.trapz(self.Ex_DDa) * np.trapz(self.Em_DAa) if np.trapz(self.Ex_DDa) * np.trapz(self.Em_DAa) else 1e-10
        self.x6 = np.trapz(self.Ex_AAd) * np.trapz(self.Em_AAd) if np.trapz(self.Ex_AAd) * np.trapz(self.Em_AAd) else 1e-10
        self.x7 = (np.trapz(self.Em_AAa) - np.trapz(self.Em_AAd)) * np.trapz(self.Ex_AAd) if (np.trapz(self.Em_AAa) - np.trapz(self.Em_AAd)) * np.trapz(self.Ex_AAd) else 1e-10
        self.x8 = np.trapz(self.Ex_AAa) * np.trapz(self.Em_AAa) if np.trapz(self.Ex_AAa) * np.trapz(self.Em_AAa) else 1e-10

        if not quiet:
            print("X0 = {:.2e}\nX1 = {:.2e}\nX2 = {:.2e}".format(self.x0,self.x1,self.x2))
            print("X3 = {:.2e}\nX4 = {:.2e}\nX5 = {:.2e}".format(self.x3,self.x4,self.x5))
            print("X6 = {:.2e}\nX7 = {:.2e}\nX8 = {:.2e}".format(self.x6,self.x7,self.x8))

        if plot:

            d_max_pos = self.d.emi_spec.argmax()
            a_max_pos = self.a.emi_spec.argmax()
            self.d_max_wl = self.wref[d_max_pos]
            self.a_max_wl = self.wref[a_max_pos]

            fig, axes = plt.subplots(3,3, figsize=(12,8))
            axes[0,0].set_title('X0 = Ex_DDd * Em_DDd = {:.2e}'.format(self.x0), pad=20)
            axes[0,0].plot(self.wref, self.Ex_DDd_graph, label='Ex_DDd')
            axes[0,0].fill_between(self.wref, self.Ex_DDd_graph, alpha=0.5)
            axes[0,0].text(self.wref[self.Ex_DDd_graph.argmax()], self.Ex_DDd_graph.max()/2, "{:.2e}".format(np.trapz(self.Ex_DDd)), horizontalalignment='center', color="black")
            axes[0,0].plot(self.wref, self.Em_DDd_graph, label='Em_DDd')
            axes[0,0].fill_between(self.wref, self.Em_DDd_graph, alpha=0.5)
            axes[0,0].text(self.wref[self.Em_DDd_graph.argmax()], self.Em_DDd_graph.max()/2, "{:.2e}".format(np.trapz(self.Em_DDd)), horizontalalignment='center', color="black")

            axes[0,1].set_title('X1 = (Em_DDa - Em_DDd) * Ex_DDd = {:.2e}'.format(self.x1), pad=20)
            axes[0,1].plot(self.wref, self.Em_DDa_graph, label='Em_DDa')
            axes[0,1].fill_between(self.wref, self.Em_DDa_graph, alpha=0.5)
            axes[0,1].text(self.wref[self.Em_DDa_graph.argmax()], self.Em_DDa_graph.max()/2, "{:.2e}".format(np.trapz(self.Em_DDa)), horizontalalignment='center', color="black")
            axes[0,1].plot(self.wref, self.Em_DDd_graph, label='Em_DDd')
            axes[0,1].fill_between(self.wref, self.Em_DDd_graph, alpha=0.5)
            axes[0,1].text(self.wref[self.Em_DDd_graph.argmax()], self.Em_DDd_graph.max()/2, "{:.2e}".format(np.trapz(self.Em_DDd)), horizontalalignment='center', color="black")
            axes[0,1].plot(self.wref, self.Ex_DDd_graph, label='Ex_DDd')
            axes[0,1].fill_between(self.wref, self.Ex_DDd_graph, alpha=0.5)
            axes[0,1].text(self.wref[self.Ex_DDd_graph.argmax()], self.Ex_DDd_graph.max()/2, "{:.2e}".format(np.trapz(self.Ex_DDd)), horizontalalignment='center', color="black")            

            axes[0,2].set_title('X2 = Ex_DDa * Em_DDa = {:.2e}'.format(self.x2), pad=20)
            axes[0,2].plot(self.wref, self.Ex_DDa_graph, label='Ex_DDa')
            axes[0,2].fill_between(self.wref, self.Ex_DDa_graph, alpha=0.5)
            axes[0,2].text(self.wref[self.Ex_DDa_graph.argmax()], self.Ex_DDa_graph.max()/2, "{:.2e}".format(np.trapz(self.Ex_DDa)), horizontalalignment='center', color="black")
            axes[0,2].plot(self.wref, self.Em_DDa_graph, label='Em_DDa')
            axes[0,2].fill_between(self.wref, self.Em_DDa_graph, alpha=0.5)
            axes[0,2].text(self.wref[self.Em_DDa_graph.argmax()], self.Em_DDa_graph.max()/2, "{:.2e}".format(np.trapz(self.Em_DDa)), horizontalalignment='center', color="black")            

            axes[1,0].set_title('X3 = Ex_DDd * Em_DAd = {:.2e}'.format(self.x3), pad=20)
            axes[1,0].plot(self.wref, self.Ex_DDd_graph, label='Ex_DDd')
            axes[1,0].fill_between(self.wref, self.Ex_DDd_graph, alpha=0.5)
            axes[1,0].text(self.wref[self.Ex_DDd_graph.argmax()], self.Ex_DDd_graph.max()/2, "{:.2e}".format(np.trapz(self.Ex_DDd)), horizontalalignment='center', color="black")
            axes[1,0].plot(self.wref, self.Em_DAd_graph, label='Em_DAd')
            axes[1,0].fill_between(self.wref, self.Em_DAd_graph, alpha=0.5)
            axes[1,0].text(self.wref[self.Em_DAd_graph.argmax()], self.Em_DAd_graph.max()/2, "{:.2e}".format(np.trapz(self.Em_DAd)), horizontalalignment='center', color="black")            

            axes[1,1].set_title('X4 = (Em_DAa - Em_DAd) * Ex_DDd = {:.2e}'.format(self.x4), pad=20)
            axes[1,1].plot(self.wref, self.Em_DAa_graph, label='Em_DAa')
            axes[1,1].fill_between(self.wref, self.Em_DAa_graph, alpha=0.5)
            axes[1,1].text(self.wref[self.Em_DAa_graph.argmax()], self.Em_DAa_graph.max()/2, "{:.2e}".format(np.trapz(self.Em_DAa)), horizontalalignment='center', color="black")                        
            axes[1,1].plot(self.wref, self.Em_DAd_graph, label='Em_DAd')
            axes[1,1].fill_between(self.wref, self.Em_DAd_graph, alpha=0.5)
            axes[1,1].text(self.wref[self.Em_DAd_graph.argmax()], self.Em_DAd_graph.max()/2, "{:.2e}".format(np.trapz(self.Em_DAd)), horizontalalignment='center', color="black")                        
            axes[1,1].plot(self.wref, self.Ex_DDd_graph, label='Ex_DDd')
            axes[1,1].fill_between(self.wref, self.Ex_DDd_graph, alpha=0.5)
            axes[1,1].text(self.wref[self.Ex_DDd_graph.argmax()], self.Ex_DDd_graph.max()/2, "{:.2e}".format(np.trapz(self.Ex_DDd)), horizontalalignment='center', color="black")                        

            axes[1,2].set_title('X5 = Ex_DDa * Em_DAa = {:.2e}'.format(self.x5), pad=20)
            axes[1,2].plot(self.wref, self.Ex_DDa_graph, label='Ex_DDa')
            axes[1,2].fill_between(self.wref, self.Ex_DDa_graph, alpha=0.5)
            axes[1,2].text(self.wref[self.Ex_DDa_graph.argmax()], self.Ex_DDa_graph.max()/2, "{:.2e}".format(np.trapz(self.Ex_DDa)), horizontalalignment='center', color="black")
            axes[1,2].plot(self.wref, self.Em_DAa_graph, label='Em_DAa')
            axes[1,2].fill_between(self.wref, self.Em_DAa_graph, alpha=0.5)
            axes[1,2].text(self.wref[self.Em_DAa_graph.argmax()], self.Em_DAa_graph.max()/2, "{:.2e}".format(np.trapz(self.Em_DAa)), horizontalalignment='center', color="black")

            axes[2,0].set_title('X6 = Ex_AAd * Em_AAd = {:.2e}'.format(self.x6), pad=20)
            axes[2,0].plot(self.wref, self.Ex_AAd_graph, label='Ex_AAd')
            axes[2,0].fill_between(self.wref, self.Ex_AAd_graph, alpha=0.5)
            axes[2,0].text(self.wref[self.Ex_AAd_graph.argmax()], self.Ex_AAd_graph.max()/2, "{:.2e}".format(np.trapz(self.Ex_AAd)), horizontalalignment='center', color="black")
            axes[2,0].plot(self.wref, self.Em_AAd_graph, label='Em_AAd')
            axes[2,0].fill_between(self.wref, self.Em_AAd_graph, alpha=0.5)
            axes[2,0].text(self.wref[self.Em_AAd_graph.argmax()], self.Em_AAd_graph.max()/2, "{:.2e}".format(np.trapz(self.Em_AAd)), horizontalalignment='center', color="black")            

            axes[2,1].set_title('X7 = (Em_AAa - Em_AAd) * Ex_AAd = {:.2e}'.format(self.x7), pad=20)
            axes[2,1].plot(self.wref, self.Em_AAa_graph, label='Em_AAa')
            axes[2,1].fill_between(self.wref, self.Em_AAa_graph, alpha=0.5)
            axes[2,1].text(self.wref[self.Em_AAa_graph.argmax()], self.Em_AAa_graph.max()/2, "{:.2e}".format(np.trapz(self.Em_AAa)), horizontalalignment='center', color="black")            
            axes[2,1].plot(self.wref, self.Em_AAd_graph, label='Em_AAd')
            axes[2,1].fill_between(self.wref, self.Em_AAd_graph, alpha=0.5)
            axes[2,1].text(self.wref[self.Em_AAd_graph.argmax()], self.Em_AAd_graph.max()/2, "{:.2e}".format(np.trapz(self.Em_AAd)), horizontalalignment='center', color="black")                        
            axes[2,1].plot(self.wref, self.Ex_AAd_graph, label='Ex_AAd')
            axes[2,1].fill_between(self.wref, self.Ex_AAd_graph, alpha=0.5)
            axes[2,1].text(self.wref[self.Ex_AAd_graph.argmax()], self.Ex_AAd_graph.max()/2, "{:.2e}".format(np.trapz(self.Ex_AAd)), horizontalalignment='center', color="black")                                    

            axes[2,2].set_title('X8 = Ex_AAa * Em_AAa = {:.2e}'.format(self.x8), pad=20)
            axes[2,2].plot(self.wref, self.Ex_AAa_graph, label='Ex_AAa')
            axes[2,2].fill_between(self.wref, self.Ex_AAa_graph, alpha=0.5)
            axes[2,2].text(self.wref[self.Ex_AAa_graph.argmax()], self.Ex_AAa_graph.max()/2, "{:.2e}".format(np.trapz(self.Ex_AAa)), horizontalalignment='center', color="black")                                                
            axes[2,2].plot(self.wref, self.Em_AAa_graph, label='Em_AAa')
            axes[2,2].fill_between(self.wref, self.Em_AAa_graph, alpha=0.5)
            axes[2,2].text(self.wref[self.Em_AAa_graph.argmax()], self.Em_AAa_graph.max()/2, "{:.2e}".format(np.trapz(self.Em_AAa)), horizontalalignment='center', color="black")                                                

            for ax in axes.flat:
                ax.set_xlim(self.d_max_wl-100, self.a_max_wl+100)
                ax.legend()

            plt.tight_layout()

    def calc_efficiency_full(self, quiet=False, filtervalues=True, shiftvalues=True, plot=False):
        self.calc_correction_factors_full(quiet=True, plot=False)
        self.df_corr_full = pd.DataFrame() 
        self.df_corr_full['time'] = self.dd_raw['time']
        self.df_corr_full['label'] = self.dd_raw['label']
        self.df_corr_full['unique_cell'] = self.dd_raw['unique_cell']

        self.df_corr_full['d_conc'] = np.divide((self.x2 * self.x4 - self.x5 * self.x1) * (self.x2 * self.aa_raw['value'] - self.x8 * self.dd_raw['value']) - (self.x2 * self.x7 - self.x8 * self.x1) * (self.x2 * self.da_raw['value'] - self.x5 * self.dd_raw['value']),
                                        (self.x2 * self.x4 - self.x5 * self.x1) * (self.x2 * self.x6 - self.x8 * self.x0) - (self.x2 * self.x7 - self.x8 * self.x1) * (self.x2 * self.x3 - self.x5 * self.x0))

        self.df_corr_full['fret_eff'] = np.divide((self.x2 * self.da_raw['value'] - self.x5 * self.dd_raw['value']) - (self.x2 * self.x3 - self.x5 * self.x0) * self.df_corr_full['d_conc'] ,
                                        (self.x2 * self.x4 - self.x5 * self.x1) * self.df_corr_full['d_conc'])
        self.df_corr_full['a_conc'] = np.divide(self.dd_raw['value'] - self.x0 * self.df_corr_full['d_conc'] - self.x1 * self.df_corr_full['d_conc'] * self.df_corr_full['fret_eff'],
                                        self.x2)

        self.data_to_use = self.df_corr_full

        if not quiet: print("Corrected Traces were calculated")

        if filtervalues:
            df_complete = self.df_corr_full
            unique_cells = df_complete['unique_cell'].unique()
            unique_cells_to_discard = []
            unique_cells_to_keep = []
            for cell in unique_cells:
                current_df = df_complete[df_complete['unique_cell'] == cell]
                if ((current_df.fret_eff > 1) | (current_df.fret_eff < 0)).any():
                    unique_cells_to_discard.append(cell)
                else:
                    unique_cells_to_keep.append(cell)
            if not quiet: print("I kept {} unique cells and discarded {} unique cells".format(len(unique_cells_to_keep), len(unique_cells_to_discard)))
            self.df_corr_full_filter = df_complete[df_complete.unique_cell.isin(unique_cells_to_keep)]
            self.data_to_use = self.df_corr_full_filter

        if shiftvalues:
            df_first_point = self.data_to_use[self.data_to_use.time==0].reset_index().drop(columns='index')
            df_first_point['diff'] = 1 - df_first_point['fret_eff']
            df_first_point = df_first_point[['unique_cell','diff']]
            self.data_to_use = self.data_to_use.merge(df_first_point, on='unique_cell', how='left')
            self.data_to_use['fret_eff_shift'] = self.data_to_use['fret_eff'] + self.data_to_use['diff']  

        if plot:
            if shiftvalues == False:
                fig, axes = plt.subplots(1,1, figsize=(8,4))
                g1 = sns.lineplot(data=self.data_to_use, x='time', y='fret_eff', ax=axes).set(title='FRET efficiency')
            elif shiftvalues:
                fig, axes = plt.subplots(1,2, figsize=(13,4))
                g1 = sns.lineplot(data=self.data_to_use, x='time', y='fret_eff', ax=axes[0]).set(title='FRET efficiency')
                g2 = sns.lineplot(data=self.data_to_use, x='time', y='fret_eff_shift', ax=axes[1]).set(title='Shifted FRET efficiency')

            
    def calc_correction_factors(self, quiet=False, plot=False):

        self.Em_DAd = self.da_exposure * self.d_dm * self.a_emi_filter * self.d.emi_spec_norm * self.camera # how much donor is in FRET channel?
        self.Em_DDd = self.dd_exposure * self.d_dm * self.d_emi_filter * self.d.emi_spec_norm * self.camera # how much donor is in DONOR channel?

        self.alpha = np.divide(np.trapz(self.Em_DAd),
                                np.trapz(self.Em_DDd))
        if not quiet: print("Donor bleedthrough alpha = {:.5f}".format(self.alpha))

        
        DAa_all_filters =  self.d_dm_ref * self.d_led * self.d_exc_filter # how much acceptor is excited for FRET channel?
        AAa_all_filters =  self.a_dm_ref * self.a_led * self.a_exc_filter  # how much acceptor is excited for ACCEPTOR channel?

        DAa_all_filters_norm = (self.d_pow * DAa_all_filters * self.wref[DAa_all_filters.argmax()]) / np.trapz(DAa_all_filters)
        AAa_all_filters_norm = (self.a_pow * AAa_all_filters * self.wref[AAa_all_filters.argmax()]) / np.trapz(AAa_all_filters)

        self.Ex_DAa = self.da_exposure * self.a.epsilon * DAa_all_filters_norm * self.a.exc_spec 
        self.Ex_AAa = self.aa_exposure * self.a.epsilon * AAa_all_filters_norm * self.a.exc_spec

        self.delta = np.divide(np.trapz(self.Ex_DAa),
                                np.trapz(self.Ex_AAa))

        if not quiet: print("Acceptor cross-excitation delta = {:.5f}".format(self.delta))

        self.Em_DAa = self.aa_exposure * self.d_dm * self.a_emi_filter * self.a.emi_spec_norm * self.camera # how much acceptor is in FRET channel?

        self.gamma = np.divide(np.trapz(self.Em_DAa), 
                                np.trapz(self.Em_DDd))
        if not quiet: print("Detection correction gamma = {:.5f}".format(self.gamma))

        # self.Ex_DDd = self.d.epsilon * self.d_dm_ref * self.d_led * self.d_exc_filter * self.d_pow * self.d.exc_spec # how much donor is excited for DONOR channel?

        # self.beta = np.divide(np.trapz(self.Ex_AAa), 
        #                         np.trapz(self.Ex_DDd))
        # if not quiet: print("Excitation correction beta = {:.2f}".format(self.beta))

        if plot:
            fig, axes = plt.subplots(1,3, figsize=(15,4))

            # Plot relevant spectra for alpha
            axes[0].plot(self.wref, self.Em_DAd, label='Em DAd')
            axes[0].plot(self.wref, self.Em_DDd, label='Em DDd')
            axes[0].fill_between(self.wref, self.Em_DAd, alpha=0.5)
            axes[0].fill_between(self.wref, self.Em_DDd, alpha=0.5)

            axes[0].text(self.wref[self.Em_DDd.argmax()], self.Em_DDd.max()/2, "{:.2e}".format(np.trapz(self.Em_DDd)), horizontalalignment='center', color="black")
            axes[0].text(self.wref[self.Em_DAd.argmax()], self.Em_DAd.max()/2, "{:.2e}".format(np.trapz(self.Em_DAd)), horizontalalignment='center', color="black")
            axes[0].set_xlim(self.wref[self.Em_DDd.argmax()]-100, self.wref[self.Em_DAd.argmax()]+100)

            axes[0].set_title(r"Donor bleedthrough $\alpha$")

            # Plot relevant spectra for delta
            axes[1].plot(self.wref, self.Ex_DAa, label='Ex DAa')
            axes[1].plot(self.wref, self.Ex_AAa, label='Ex AAa')
            axes[1].fill_between(self.wref, self.Ex_DAa, alpha=0.5)
            axes[1].fill_between(self.wref, self.Ex_AAa, alpha=0.5)

            axes[1].text(self.wref[self.Ex_DAa.argmax()], self.Ex_DAa.max()/2, "{:.2e}".format(np.trapz(self.Ex_DAa)), horizontalalignment='center', color="black")
            axes[1].text(self.wref[self.Ex_AAa.argmax()], self.Ex_AAa.max()/2, "{:.2e}".format(np.trapz(self.Ex_AAa)), horizontalalignment='center', color="black")
            axes[1].set_xlim(self.wref[self.Ex_DAa.argmax()]-100, self.wref[self.Ex_AAa.argmax()]+100)

            axes[1].set_title(r"Acceptor cross-excitation $\delta$")

            # Plot relevant spectra for gamma
            axes[2].plot(self.wref, self.Em_DAa, label='Em DAa')
            axes[2].plot(self.wref, self.Em_DDd, label='Em DDd')
            axes[2].fill_between(self.wref, self.Em_DAa, alpha=0.5)
            axes[2].fill_between(self.wref, self.Em_DDd, alpha=0.5)

            axes[2].text(self.wref[self.Em_DDd.argmax()], self.Em_DDd.max()/2, "{:.2e}".format(np.trapz(self.Em_DDd)), horizontalalignment='center', color="black")
            axes[2].text(self.wref[self.Em_DAa.argmax()], self.Em_DAa.max()/2, "{:.2e}".format(np.trapz(self.Em_DAa)), horizontalalignment='center', color="black")
            axes[2].set_xlim(self.wref[self.Em_DDd.argmax()]-100, self.wref[self.Em_DAa.argmax()]+100)
            axes[2].set_title(r"Detection correction $\gamma$")

            # # Plot relevant spectra for beta
            # axes[3].plot(self.wref, self.Ex_DDd, label='Ex DDd')
            # axes[3].plot(self.wref, self.Ex_AAa, label='Ex AAa')
            # axes[3].fill_between(self.wref, self.Ex_DDd, alpha=0.5)
            # axes[3].fill_between(self.wref, self.Ex_AAa, alpha=0.5)

            # axes[3].text(self.wref[self.Ex_DDd.argmax()], self.Ex_DDd.max()/2, "{:.2f}".format(np.trapz(self.Ex_DDd)), horizontalalignment='center', color="black")
            # axes[3].text(self.wref[self.Ex_AAa.argmax()], self.Ex_AAa.max()/2, "{:.2f}".format(np.trapz(self.Ex_AAa)), horizontalalignment='center', color="black")
            # axes[3].set_xlim(self.wref[self.Ex_DDd.argmax()]-100, self.wref[self.Ex_AAa.argmax()]+100)
            # axes[3].set_title(r"Excitation correction $\beta$")

            for ax in axes.flat:
                ax.legend()

            plt.tight_layout()

    def calc_efficiency(self, quiet=False, filtervalues=True, shiftvalues=True, plot=False):
        # self.calc_correction_factors(quiet=True, plot=False)
        self.df_corrections = pd.DataFrame() 
        self.df_corrections['time'] = self.dd_raw['time']
        self.df_corrections['label'] = self.dd_raw['label']
        self.df_corrections['unique_cell'] = self.dd_raw['unique_cell']
        self.df_corrections['fc'] = self.da_raw['value'] - (self.alpha * self.dd_raw['value']) - (self.delta * self.aa_raw['value'])
        self.df_corrections['fret_eff'] = 1 / (np.divide((self.gamma * self.dd_raw['value']), self.df_corrections['fc']) + 1)
    
        self.data_to_use = self.df_corrections

        if not quiet: print("Efficiency was calculated")

        if filtervalues:
            df_complete = self.df_corrections
            unique_cells = df_complete['unique_cell'].unique()
            unique_cells_to_discard = []
            unique_cells_to_keep = []
            for cell in unique_cells:
                current_df = df_complete[df_complete['unique_cell'] == cell]
                if ((current_df.fret_eff > 1) | (current_df.fret_eff < 0)).any():
                    unique_cells_to_discard.append(cell)
                else:
                    unique_cells_to_keep.append(cell)
            if not quiet: print("I kept {} unique cells and discarded {} unique cells".format(len(unique_cells_to_keep), len(unique_cells_to_discard)))
            self.df_corrections_filter = df_complete[df_complete.unique_cell.isin(unique_cells_to_keep)]
            self.data_to_use = self.df_corrections_filter

        if shiftvalues:
            df_first_point = self.data_to_use[self.data_to_use.time==0].reset_index().drop(columns='index')
            df_first_point['diff'] = 1 - df_first_point['fret_eff']
            df_first_point = df_first_point[['unique_cell','diff']]
            self.data_to_use = self.data_to_use.merge(df_first_point, on='unique_cell', how='left')
            self.data_to_use['fret_eff_shift'] = self.data_to_use['fret_eff'] + self.data_to_use['diff'] 

        if plot:
            if not shiftvalues:
                fig, axes = plt.subplots(1,2, figsize=(15,4))
                g1 = sns.lineplot(data=self.data_to_use, x='time', y='fc', ax=axes[0]).set(title='Sensitized emission')
                g2 = sns.lineplot(data=self.data_to_use, x='time', y='fret_eff', ax=axes[1]).set(title='FRET efficiency')
            else:
                fig, axes = plt.subplots(1,3, figsize=(15,4))
                g1 = sns.lineplot(data=self.data_to_use, x='time', y='fc', ax=axes[0]).set(title='Sensitized emission')
                g2 = sns.lineplot(data=self.data_to_use, x='time', y='fret_eff', ax=axes[1]).set(title='FRET efficiency')
                g3 = sns.lineplot(data=self.data_to_use, x='time', y='fret_eff_shift', ax=axes[2]).set(title='Shifted FRET efficiency')

    # def calc_stoichiometry(self, plot=False):
    #     try:
    #        self.df_corrections['stoichiometry'] = np.divide(self.gamma*self.dd_raw['value'] + self.da_raw['value'],self.gamma*self.dd_raw['value']+self.da_raw['value']+(self.aa_raw['value']/self.beta))
    #        print("Stoichiometry was calculated")
    #     except NameError:
    #         self.df_corrections = pd.DataFrame() 
    #         self.df_corrections['time'] = self.dd_raw['time']
    #         self.df_corrections['label'] = self.dd_raw['label']
    #         self.df_corrections['stoichiometry'] = np.divide(self.gamma*self.dd_raw['value'] + self.da_raw['value'],
    #         self.gamma*self.dd_raw['value']+self.da_raw['value']+ (self.aa_raw['value']/self.beta))
    #         print("Stoichiometry was calculated")

    #     if plot:
    #         fig, axes = plt.subplots(figsize=(10,6))
    #         g1 = sns.lineplot(data=self.df_corrections, x='time', y='stoichiometry', ax=axes).set(title='Stoichiometry')

    def plot_raw_data(self):    
        fig, axes = plt.subplots(1,3, figsize=(15,4))
        g1 = sns.lineplot(data=self.dd_raw, x='time', y='value', ax=axes[0]).set(title='Donor RAW')
        g2 = sns.lineplot(data=self.da_raw, x='time', y='value', ax=axes[1]).set(title='FRET RAW')
        g3 = sns.lineplot(data=self.aa_raw, x='time', y='value', ax=axes[2]).set(title='Acceptor RAW')

    def compare_corrections(self, filtervalues=True, shift=True):
        if shift:
            traces_to_plot = 'fret_eff_shift'
        else: traces_to_plot = 'fret_eff'
        fig, axes = plt.subplots(figsize=(6,4))
        self.calc_correction_factors(quiet=True)
        self.calc_efficiency(quiet=True, filtervalues=filtervalues, shiftvalues=True, plot=False)
        sns.lineplot(data=self.data_to_use, x='time', y=traces_to_plot, ax=axes, label='Traditional Corrections')
        self.calc_correction_factors_full(quiet=True)
        self.calc_efficiency_full(quiet=True, filtervalues=filtervalues, shiftvalues=True, plot=False)
        sns.lineplot(data=self.data_to_use, x='time', y=traces_to_plot, ax=axes, label='Complete Corrections')
        plt.title('Corrections comparison for ' + self.name)
        if not shift: plt.ylim(0,1)
        plt.legend()
        plt.tight_layout()