"""
The following code is a modified version of scFEA 
(https://github.com/changwn/scFEAAlghamdi et al, Genome Research, 2021) 
which includes exchange reactions between cell types.
"""
from torch.utils.data import Dataset

class MyDataset(Dataset):
    def __init__(self, data, label, info, transform = None):
        self.data = data
        self.label = label
        self.info = info
        self.transform = transform

    def __getitem__(self, index):
        x = self.data[index]
        y = self.label[index]
        z = self.info[index]
        if self.transform:
            x = self.transform(x)
        return x, y, z

    def __len__(self):
        return len(self.data)