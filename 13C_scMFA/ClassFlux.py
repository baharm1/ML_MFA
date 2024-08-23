import torch
import torch.nn as nn


class FLUX(nn.Module):
    def __init__(self, n_modules, n_genes):
        super(FLUX, self).__init__()
        # gene to flux
        self.n_genes = n_genes
        self.n_modules = n_modules
        
        self.gene_features = nn.ModuleList([nn.Sequential(
            nn.Linear(self.n_genes, 8, bias=False),
            nn.Tanhshrink(),
            nn.Linear(8, 1, bias=False)) for i in range(self.n_modules)])
                

    def forward(self, x, module_scale, cmMat, mmMat):

        for i in range(self.n_modules):
            
            x_block = x[:, i * self.n_genes : (i + 1) * self.n_genes]
            gene_subnet = self.gene_features[i]
            
            if i == 0:
                gene_flux = gene_subnet(x_block)
            else:
                gene_flux = torch.cat((gene_flux, gene_subnet(x_block)), 1)

        compound_acc_gene = gene_flux @ cmMat.T
        
        penalty_module_scale = 1 / (0.01 + module_scale)
        penalized_mid_bal = (penalty_module_scale * gene_flux) @ mmMat.T
        
        isotopologue_acc = gene_flux @ mmMat.T #matmul
        
        
        return gene_flux, compound_acc_gene, isotopologue_acc, penalized_mid_bal