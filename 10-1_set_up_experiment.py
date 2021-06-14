# -*- coding: utf-8 -*-
"""
Created on Tue Nov  3 20:02:25 2020

@author: user
"""
import TisOpt

#%% Optimize 1
opt = TisOpt.TissueOptimizer("Lung", n_pool=1000)
# https://www.uniprot.org/uniprot/C5MKY7
egfp = "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK"
opt.optimize(egfp)
best_egfp_lung = opt.select_best(by={"MFE":"min","MFEini":"max","CAI":"max","CPB":"max","ENC":"min"},homopolymers=7, top=10)
# https://www.uniprot.org/uniprot/X5DSL3
mcherry = "MVSKGEEDNMAIIKEFMRFKVHMEGSVNGHEFEIEGEGEGRPYEGTQTAKLKVTKGGPLPFAWDILSPQFMYGSKAYVKHPADIPDYLKLSFPEGFKWERVMNFEDGGVVTVTQDSSLQDGEFIYKVKLRGTNFPSDGPVMQKKTMGWEASSERMYPEDGALKGEIKQRLKLKDGGHYDAEVKTTYKAKKPVQLPGAYNVNIKLDITSHNEDYTIVEQYERAEGRHSTGGMDELYK"
opt.optimize(mcherry)
best_mcherry_lung = opt.select_best(by={"MFE":"min","MFEini":"max","CAI":"max","CPB":"max","ENC":"min"},homopolymers=7, top=10)

#%% Optimize 2
opt = TisOpt.TissueOptimizer("Kidney", n_pool=1000)
# https://www.uniprot.org/uniprot/C5MKY7
egfp = "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK"
opt.optimize(egfp)
best_egfp_kidney = opt.select_best(by={"MFE":"min","MFEini":"max","CAI":"max","CPB":"max","ENC":"min"},homopolymers=7, top=10)
# https://www.uniprot.org/uniprot/X5DSL3
mcherry = "MVSKGEEDNMAIIKEFMRFKVHMEGSVNGHEFEIEGEGEGRPYEGTQTAKLKVTKGGPLPFAWDILSPQFMYGSKAYVKHPADIPDYLKLSFPEGFKWERVMNFEDGGVVTVTQDSSLQDGEFIYKVKLRGTNFPSDGPVMQKKTMGWEASSERMYPEDGALKGEIKQRLKLKDGGHYDAEVKTTYKAKKPVQLPGAYNVNIKLDITSHNEDYTIVEQYERAEGRHSTGGMDELYK"
opt.optimize(mcherry)
best_mcherry_kidney = opt.select_best(by={"MFE":"min","MFEini":"max","CAI":"max","CPB":"max","ENC":"min"},homopolymers=7, top=10)

#%% Write
#best_egfp_lung.to_csv("../results/egfp_lung_optimized.csv")
#best_egfp_kidney.to_csv("../results/egfp_kidney_optimized.csv")
#best_mcherry_lung.to_csv("../results/mcherry_lung_optimized.csv")
#best_mcherry_kidney.to_csv("../results/mcherry_kidney_optimized.csv")