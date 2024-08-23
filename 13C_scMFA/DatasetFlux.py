# -*- coding: utf-8 -*-
"""
@author: wnchang pulled y1zhou
"""

from torch.utils.data import Dataset


class MyDataset(Dataset):
    def __init__(self, geneExpr, geneExprScale, module_scale, transform=None):
        self.geneExpr = geneExpr
        self.geneExprScale = geneExprScale
        self.module_scale = module_scale
        
        self.transform = transform

    def __getitem__(self, index):
        x = self.geneExpr[index]
        y = self.geneExprScale[index]
        z = self.module_scale[index]
        
        if self.transform:
            x = self.transform(x)
        return x, y, z

    def __len__(self):
        return len(self.geneExpr)