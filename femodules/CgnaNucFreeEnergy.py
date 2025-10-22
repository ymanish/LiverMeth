####This module calculates nucleosome free energies using the CGNA Parameterization using the non-optimized code.
####This also used the marganalization of the Cgna Parameters to 6 degrees of freedom per base pair step, similar to RBP model.
####This modules is temporary and will be replaced by optimized code in the future, where are the DOF are included.
import os
# #Only import env_settings if the flag is set (default is to import)
if os.environ.get("IMPORT_ENV_SETTINGS", "1") == "1":
    from src.config.env_settings import *  # Triggers env_settings import

# from py_analysis.config.env_settings import *  # Triggers env_settings import
import numpy as np # Inherits the thread limits
from typing import Optional

from src.config.custom_types import FreeEnergyResult
from femodules import NUC_STATE_PATH, K_POSRESC_PATH


from methods import nucleosome_free_energy, read_nucleosome_triads, GenStiffness, calculate_midstep_triads
from binding_model import binding_model_free_energy
# from line_profiler import profile 
# from src.core.helper.fn_wraps import timeit


class NucleosomeBreath:

    def __init__(self,  
                nuc_method:str='crystal', free_dna_method: Optional[str] = None):
        
        self.nuc_method = nuc_method
        self.free_dna_method = free_dna_method

        print(f"Using nucleosome method: {self.nuc_method}")
        print(f"Using free DNA method: {self.free_dna_method }")

        self.genstiff_nuc = GenStiffness(method=self.nuc_method )

        if self.free_dna_method  is not None:
            ### Used in the case of multiharmonic parameterization
            self.genstiff_freedna = GenStiffness(method=self.free_dna_method )

        self.triadfn = NUC_STATE_PATH
        self.nuctriads = read_nucleosome_triads(self.triadfn)

        self.fn = K_POSRESC_PATH
        self.K_resc = np.load(self.fn)

        self.nuc_mu0 = calculate_midstep_triads(triad_ids = self._select_phosphate_bind_sites(), 
                                           nucleosome_triads = self.nuctriads)



    def _select_phosphate_bind_sites(self, left=0, right=13):

        phosphate_bind_sites = [2, 6, 14, 17, 24, 29, 34, 38, 
                                    45, 49, 55, 59, 65, 69, 76, 
                                    80, 86, 90, 96, 100, 107, 111, 
                                    116, 121, 128, 131, 139, 143]
        
        return phosphate_bind_sites[left*2:(right*2)+2]

    def _get_left_right_open(self, left:int, right:int, style:str="b_index") -> tuple:
        if style == "b_index":
            ## Inputs are bound‐site indices (0-13) 
            ## Each index covers two phosphate positions:
            ## bound site i <-> phosphates (2*i, 2*i+1). 
            ## this converts the bound index to how many open phosphate sites required to be there
            ## Example: left=1, right=11, means l_open=2, r_open=4 so you need 2 open phosphate sites on left and 4 on right
            l_open = 2*left
            r_open  = 28-(2*right)-2

        elif style == "ph_index":
            ## Inputs are phosphate‐site indices directly (0-27).
            ## Same as the bi style but directly at the phosphate sites level
            ## Example: left=1, right=26, means l_open=1, r_open=1
            l_open = left
            r_open  = 28-(right)-1
           
        elif style == "open_sites":
            ## Inputs already specify exactly how many phosphates are open
            ## on each side—no conversion needed.
            ## Example: left=0, right=27, means l_open=0, r_open=27
            l_open = left
            r_open  = right

        else:
            raise ValueError("Invalid style. Use 'b_index' or 'open_sites' or 'ph_index'.")
        
        ### The left and right open are the number of open phosphates on each side
        return l_open, r_open

    def calculate_free_energy_soft(self, seq601:str, left:int, right:int, 
                                   id:Optional[str]=None, subid:Optional[str]=None,
                                    kresc_factor:float = 1, style:str="b_index", bound_ends:str="exclude")-> FreeEnergyResult:
        
        """ bounds_ends: str = 'include' or 'exclude' , this is used when the dna is full bound on the histone core, 
                                in that case there two outermost base pair steps which can either be calculated as 
                                bound or unbound separately. Or you can include them as just bound dna part of the histone core"""
        
        stiff, gs = self.genstiff_nuc.gen_params(seq601, use_group=True, sparse=False)

        l_open, r_open = self._get_left_right_open(left, right, style)


        F_dict = binding_model_free_energy(
            gs,
            stiff,    
            self.nuc_mu0,
            self.K_resc*kresc_factor,
            left_open=l_open,
            right_open=r_open,
            use_correction=True,
        )

        # Kentries = [1,1,1,10,10,10]
        # diags = np.concatenate([Kentries]*len(nuc_mu0))
        # K = np.diag(diags)

        # F_dict = binding_model_free_energy(
        #     gs,
        #     stiff,    
        #     nuc_mu0,
        #     K,
        #     left_open=l_open,
        #     right_open=r_open,
        #     use_correction=True,
        # )


        if self.free_dna_method is not None:
            bound_freedna_fe, unbound_freedna_fe, new_freedna_fe = self.calculate_bound_nonbound_dna_energy(seq=seq601, 
                                                                                      left_open=l_open, 
                                                                                      right_open=r_open,
                                                                                        bound_ends=bound_ends)
            

            F_601_new = (F_dict['F'] - F_dict['F_freedna']) +  bound_freedna_fe + unbound_freedna_fe
            F_free_new = new_freedna_fe
            F_entropy_new = F_601_new - F_dict['F_enthalpy'] 
            F_entalap = F_dict['F_enthalpy']
            
            return FreeEnergyResult(F_601_new, F_entropy_new, F_entalap, F_free_new, id, subid)



        F601 = F_dict['F']
        F_entrop = F_dict['F_entropy']
        F_entalap = F_dict['F_enthalpy']
        F_free = F_dict['F_freedna']

        return FreeEnergyResult(F601, F_entrop, F_entalap, F_free, id, subid)
    
    def calculate_bound_nonbound_dna_energy(self, seq:str, left_open:int, right_open:int, bound_ends: str = "include")-> tuple:
        """ Calculate the free energy of the bound and unbound DNA parts of the nucleosome.
        Used in the case of multiharmonic parameterization.

         bounds_ends: str = 'include' or 'exclude' , this is used when the dna is full bound on the histone core, 
                                in that case there two outermost base pair steps which can either be calculated as 
                                bound or unbound separately. Or you can include them as just bound dna part of the histone core"""

        if bound_ends not in {"include", "exclude"}:
            raise ValueError("bound_ends must be either 'include' or 'exclude'.")

        full_stiff_dna_unbound, gs_dna = self.genstiff_freedna.gen_params(seq, use_group=True)
        full_stiff_dna_bound, gs_dna = self.genstiff_nuc.gen_params(seq, use_group=True)

        if full_stiff_dna_unbound.shape[0] != full_stiff_dna_bound.shape[0]:
            raise ValueError("Stiffness matrices for unbound and bound DNA must have the same size.")

        phosphate_bind_sites = self._select_phosphate_bind_sites()

        logdet_sign, logdet = np.linalg.slogdet(full_stiff_dna_unbound)
        new_free_dna_fe = -0.5*len(full_stiff_dna_unbound)*np.log(2*np.pi) + 0.5*logdet



        if bound_ends == "include" and left_open == 0 and right_open == 0:
            logdet_sign_b, logdet_b = np.linalg.slogdet(full_stiff_dna_bound)
            Fe_bound_dna = -0.5*len(full_stiff_dna_bound)*np.log(2*np.pi) + 0.5*logdet_b
            return Fe_bound_dna, 0.0, new_free_dna_fe
        
        bound_locs = phosphate_bind_sites[left_open:len(phosphate_bind_sites)-right_open]
        
        start = 6 * bound_locs[0]
        end = 6 * (bound_locs[-1] + 1)
        bound_block_stiff = full_stiff_dna_bound[start:end, start:end]
        
        unbound_indices = np.r_[0:start, end:full_stiff_dna_bound.shape[0]] 
        unbound_stiffness = full_stiff_dna_unbound[unbound_indices][:, unbound_indices]


        logdet_sign_b, logdet_b = np.linalg.slogdet(bound_block_stiff)
        Fe_bound_dna = -0.5*len(bound_block_stiff)*np.log(2*np.pi) + 0.5*logdet_b
        
        
        logdet_sign_ub, logdet_ub = np.linalg.slogdet(unbound_stiffness)
        Fe_unbound_dna = -0.5*len(unbound_stiffness)*np.log(2*np.pi) + 0.5*logdet_ub

        return Fe_bound_dna, Fe_unbound_dna, new_free_dna_fe

    def calculate_free_energy_hard(self, seq147:str, left:int, right:int, id:Optional[str]=None, subid:Optional[str]=None)-> FreeEnergyResult:
        stiff, gs = self.genstiff_nuc.gen_params(seq147, use_group=True)



        mid_array = self._select_phosphate_bind_sites(left, right)
        F_dict = nucleosome_free_energy(gs, stiff, mid_array, self.nuctriads, use_correction=True)


        F601 = F_dict['F']
        F_entrop = F_dict['F_entropy']
        F_entalap = F_dict['F_enthalpy']
        F_free = F_dict['F_freedna']



        # return F601, F_entrop, F_entalap, F_free, F_Diff
        return FreeEnergyResult(F601, F_entrop, F_entalap, F_free, id, subid)



if __name__ == "__main__":
    import time
    start = time.perf_counter()
    Seq601 = "CTGGAGAATCCCGGTGCCGAGGCCGCTCAATTGGTCGTAGACAGCTCTAGCACCGCTTAAACGCACGTACGCGCTGTCCCCCGCGTTTTAACCGCCAAGGGGATTACTCCCTAGTCTCCAGGCACGTGTCAGATATATACATCCTGT"
    # Example usage
    nuc_breath = NucleosomeBreath(nuc_method='hybrid', free_dna_method=None)
    result = nuc_breath.calculate_free_energy_soft(seq601=Seq601, left=0, right=13, style="b_index")
    print(f"Free energy: {result.F}, Entropy: {result.F_entropy}, Enthalpy: {result.F_enthalpy}, Free DNA energy: {result.F_freedna}, dF: {result.F - result.F_freedna}")


    # free_dna_energy.get_free_dna_energy(fullseq=Seq601, left=0, right=13, style="b_index")
    # print(f"Free

    # for i in range(10):
    #     nuc_breath = NucleosomeBreath(nuc_method='hybrid', free_dna_method=None)
    #     result = nuc_breath.calculate_free_energy_soft(seq601=Seq601, left=0, right=13, style="b_index")
    #     # print(f"Left: {i}, Right: 13, Result: {result}")
    #     print(f"Free energy: {result.F}, Entropy: {result.F_entropy}, Enthalpy: {result.F_enthalpy}, Free DNA energy: {result.F_freedna}, dF: {result.F - result.F_freedna}")

    # # # for i in range(0, 14):
    #     result = nuc_breath.calculate_free_energy_soft(seq601=Seq601, left=0, right=i, style="b_index")
    #     print(f"Left: 0, Right: {i}, Result: {result}")
    

    # for i in range(0, 14):
    #     for j in range(0, 14):
    #         if j >= i:
    #             result = nuc_breath.calculate_free_energy_soft(seq601=Seq601, left=i, right=j, style="b_index")
                # print(f"L:{i}, R:{j} ------ Free energy: {result.F}, Entropy: {result.F_entropy}, Enthalpy: {result.F_enthalpy}, Free DNA energy: {result.F_freedna}, dF: {result.F - result.F_freedna}")
    # result = nuc_breath.calculate_free_energy_soft(seq601=Seq601, left=0, right=13, style="b_index")
    end = time.perf_counter()
    print(f"Time taken: {end - start:.2f} seconds")




    
    # nuc_breath.calculate_nonbound_dna_energy(seq=Seq601, left_open=0, right_open=0, bound_ends="exclude")

