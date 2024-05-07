import torch
from torch import nn
import torch.optim as optim
from torch.autograd import Variable
from torch.optim.lr_scheduler import StepLR, MultiStepLR
from sklearn.metrics import mean_squared_error
from torch.utils.data import DataLoader, TensorDataset
import pandas as pd
import numpy as np
import math

df = pd.read_csv('North_America_50kMutations_5perc_missing.csv')

# Load dataset
dna_sequences = pd.read_csv('North_America_50kMutations.tsv', delimiter='\t').iloc[:,2:].T.values
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
dna_sequences = torch.tensor(dna_sequences, dtype=torch.float32).to(device)

test_dna_sequences = pd.read_csv('North_America_50kMutations_5perc_missing.csv').iloc[:,3:].T
mask =  np.isnan(test_dna_sequences).values
test_dna_sequences = test_dna_sequences.fillna(0).values
test_dna_sequences = torch.tensor(test_dna_sequences, dtype=torch.float32).to(device)

X_train = torch.tensor(dna_sequences, dtype=torch.float).to(device)
X_test = torch.tensor(test_dna_sequences, dtype=torch.float)

train_loader = DataLoader(TensorDataset(X_train), batch_size=32, shuffle=True)

class Generator(nn.Module):
    def __init__(self, input_dim, hidden_dim):
        super(Generator, self).__init__()
        self.model = nn.Sequential(
            nn.Linear(input_dim, hidden_dim*2),
            nn.ReLU(),
            nn.Linear(hidden_dim*2, hidden_dim*2),
            nn.ReLU(),
            nn.Linear(hidden_dim*2, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, input_dim),
            nn.Sigmoid()
        )
    def forward(self, x):
        return self.model(x)

class Discriminator(nn.Module):
    def __init__(self, input_dim, hidden_dim):
        super(Discriminator, self).__init__()
        self.model = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, 1),
            nn.Sigmoid()
        )
    def forward(self, x):
        return self.model(x)

input_dim = X_train.shape[1]
hidden_dim = 256
epochs = 1000

generator = Generator(input_dim, hidden_dim).to(device)
discriminator = Discriminator(input_dim, hidden_dim).to(device)

optimizer_G = torch.optim.Adam(generator.parameters(), lr=0.0001)
optimizer_D = torch.optim.Adam(discriminator.parameters(), lr=0.0001)

criterion = nn.BCELoss()

for epoch in range(epochs): 
    #  Train Generator
    for _ in range(10):
        optimizer_G.zero_grad()
        gen_samples = generator(X_train)

        g_loss = criterion(discriminator(gen_samples), torch.ones((X_train.size(0), 1)).to(device))
        g_loss.backward()
        optimizer_G.step()
        
    #  Train Discriminator
    optimizer_D.zero_grad()

    real_loss = criterion(discriminator(X_train), torch.ones((X_train.size(0), 1)).to(device))
    fake_loss = criterion(discriminator(gen_samples.detach()), torch.zeros((X_train.size(0), 1)).to(device))
    d_loss = (real_loss + fake_loss) / 2
    d_loss.backward()
    optimizer_D.step()

    print(f'Epoch {epoch+1}/{epochs}, D-Loss: {d_loss.item():.6f}, G-Loss: {g_loss.item():.6f}')

generator.eval()
with torch.no_grad():
    X_test_imputed = generator(X_test)
    
rmse = np.sqrt(mean_squared_error(X_train.cpu(), X_test_imputed.cpu()))
print(f'RMSE on full data: {rmse:.6f}')

X_test[mask] = X_test_imputed[mask]

missing_rmse = np.sqrt(mean_squared_error(X_train.cpu(), X_test.cpu()))
print(f'RMSE on missing data: {missing_rmse:.6f}')

X_test = np.transpose(X_test.cpu(), (1,0))

df.iloc[:, 3:] = X_test.numpy()
df.to_csv('North_America_50kMutations_5perc_missing_ganfilled.csv', index=False)