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

class VAE(nn.Module):
    def __init__(self, input_dim, hidden_dim, latent_dim):
        super(VAE, self).__init__()
        self.encoder = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, latent_dim * 2)  # mean and variance to sample
        )
        self.decoder = nn.Sequential(
            nn.Linear(latent_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, input_dim),
            nn.Sigmoid() 
        )
    def reparameterize(self, mu, logvar):
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        return mu + eps * std
    def forward(self, x):
        h = self.encoder(x)
        mu, logvar = torch.chunk(h, 2, dim=1)  # split into mean and variance
        z = self.reparameterize(mu, logvar)
        return self.decoder(z), mu, logvar

def vae_loss_function(recon_x, x, mu, logvar):
    recon_loss = nn.functional.binary_cross_entropy(recon_x, x, reduction='sum')
    kld = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp())
    return (recon_loss + kld).mean()
    
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

#Hyperparameters
num_features = X_train.shape[1]
hidden_dim = 128
latent_dim = 256
num_layers = 1
epochs = 10000
model = VAE(num_features, hidden_dim, latent_dim).to(device)
optimizer = optim.AdamW(model.parameters(), lr=0.001)
scheduler = MultiStepLR(optimizer, milestones=[5000], gamma=0.1)

for epoch in range(epochs):
    train_loss = 0
    model.train()
    for batch_idx, (data,) in enumerate(train_loader):
        data = data.to(device)
        optimizer.zero_grad()
        recon_batch, mu, logvar = model(data)
        loss = vae_loss_function(recon_batch, data, mu, logvar)
        loss.backward()
        train_loss += loss.item()
        optimizer.step()
    scheduler.step()
    print(f"Epoch [{epoch+1}/{epochs}] - Loss: {math.sqrt(train_loss / len(train_loader.dataset)):.7f}")
 
# Impute missing values
model.eval()
with torch.no_grad():
    X_test_imputed, _, _ = model(X_test)        
    
rmse = np.sqrt(mean_squared_error(X_train.cpu(), X_test_imputed.cpu()))
print(f'RMSE on full data: {rmse:.6f}')

X_test[mask] = X_test_imputed[mask]

missing_rmse = np.sqrt(mean_squared_error(X_train.cpu(), X_test.cpu()))
print(f'RMSE on missing data: {missing_rmse:.6f}')

X_test = np.transpose(X_test.cpu(), (1,0))

df.iloc[:, 3:] = X_test.numpy()
df.to_csv('North_America_50kMutations_5perc_missing_vaefilled.csv', index=False)