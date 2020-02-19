from so_sz_y import SOSZPowerSpectrumLikelihood

so_sz_y = SOSZPowerSpectrumLikelihood()

info = {'likelihood': {'so_sz_y': so_sz_y},

         'params': {#'A_s': 2.215e-9,

                    'omega_cdm': {'latex': 'omega_c', 'prior': {'max': 0.15, 'min': 0.09}},

 # 'A_CIB':  # sampled! with uniform prior on [-1,1] and ref pdf N(mu=0,sigma=0.25)
#    {'prior':
#      {'min': -1,
#      'max':  1},
#    'ref':
#      {'dist': 'norm',
#      'loc':   0,
#      'scale': 0.25},
#    'latex': '\mathcal{B}_2',
#    'proposal': 0.25}
#                     'tau_reio': 0.055,
                   },
         'sampler': {'mcmc': {'covmat': 'auto', 'drag': True, 'proposal_scale': 1.9}},
         'theory': {'classy':
                       {'extra_args':
                         {
                        'omega_b' : 0.0245,
                        'A_s' : 1.884e-9,
                        'h' : 0.7,
                        'n_s' : 0.95,
                        'tau_reio' : 0.0544,

                        'HSEbias' : 1.41,

                        'output' : 'tSZ',
                        'units for tSZ spectrum' : 'dimensionless',
                        'path_to_class' : '/Users/boris/Work/CLASS-SZ/SO-SZ/class_sz',
                        'mass function' : 'M500',
                        'pressure profile' : 'A10',

                        'multipoles_sz' : 'P15',
                        'nlSZ' : 18,

                        'M1SZ' : 5e13,
                        'M2SZ' : 1e16,

                        'z1SZ' : 1e-5,
                        'z2SZ' : 4.,
                        'z_max_pk' : 4.,

                        'N_ur' : 0.00641,
                        'N_ncdm' : 1,
                        'deg_ncdm' : 3,
                        'm_ncdm' : 0.02,
                        'T_ncdm' : 0.71611,


                        'YHe' : 'BBN',
                        'k_pivot' : 0.05,
                        'sz_verbose' : 1
                            }}
                   }
       }


from cobaya.model import get_model
model = get_model(info)
print(model.loglike({}, cached=False))
