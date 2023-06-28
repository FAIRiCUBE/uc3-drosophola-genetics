import numpy as np
import math
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
import pandas as pd
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
import torch.nn.functional as F

# Load dataset
dna_sequences = pd.read_csv('North_America_50kMutations.tsv', delimiter='\t').iloc[:,2:].T.values

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

dna_sequences = torch.tensor(dna_sequences, dtype=torch.float32).to(device)
train_loader = DataLoader(TensorDataset(dna_sequences), batch_size=32, shuffle=True)

class Encoder(nn.Module):
    def __init__(self, input_size, hidden_size, latent_size):
        super(Encoder, self).__init__()
        self.fc1 = nn.Linear(input_size, hidden_size)
        self.fc_mu = nn.Linear(hidden_size, latent_size)
        self.fc_logvar = nn.Linear(hidden_size, latent_size)

    def forward(self, x):
        h = torch.relu(self.fc1(x))
        mu = self.fc_mu(h)
        logvar = self.fc_logvar(h)
        return mu, logvar

class Decoder(nn.Module):
    def __init__(self, latent_size, hidden_size, output_size):
        super(Decoder, self).__init__()
        self.fc1 = nn.Linear(latent_size, hidden_size)
        self.fc2 = nn.Linear(hidden_size, output_size)

    def forward(self, x):
        h = torch.relu(self.fc1(x))
        out = torch.sigmoid(self.fc2(h))
        return out

class VAE(nn.Module):
    def __init__(self, input_size, hidden_size, latent_size):
        super(VAE, self).__init__()
        self.encoder = Encoder(input_size, hidden_size, latent_size)
        self.decoder = Decoder(latent_size, hidden_size, input_size)

    def reparameterize(self, mu, logvar):
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        return mu + eps * std

    def forward(self, x):
        mu, logvar = self.encoder(x)
        z = self.reparameterize(mu, logvar)
        return self.decoder(z), mu, logvar

def vae_loss(recon_x, x, mu, logvar, beta=1):
    recon_loss = nn.functional.mse_loss(recon_x, x)
    kld_loss = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp())
    return recon_loss + beta * kld_loss

# Hyperparameters
input_size = 50024
hidden_size = 1024
latent_size = 512
epochs = 500
batch_size = 32
lr = 0.0001

# Initialize the VAE model and optimizer
vae = VAE(input_size, hidden_size, latent_size).to(device)
optimizer = optim.Adam(vae.parameters(), lr=lr)

# Training loop
for epoch in range(epochs):
    train_loss = 0
    for batch_idx, (data,) in enumerate(train_loader):
        data = data.to(device)
        optimizer.zero_grad()
        recon_batch, mu, logvar = vae(data)
        loss = vae_loss(recon_batch, data, mu, logvar)
        loss.backward()
        train_loss += loss.item()
        optimizer.step()

    print(f"Epoch [{epoch+1}/{epochs}] - Loss: {math.sqrt(train_loss / len(train_loader.dataset)):.4f}")
    
#Impute Missing values
test_dna_sequences = pd.read_csv('North_America_50kMutations_5perc_missing.csv').iloc[:,3:].T.fillna(0).values
test_dna_sequences = torch.tensor(test_dna_sequences, dtype=torch.float32).to(device)
test_loader = DataLoader(TensorDataset(test_dna_sequences, dna_sequences), batch_size=32, shuffle=False)

def calculate_rmse(recon_data, data):
    mse = nn.functional.mse_loss(recon_data, data, reduction='mean')
    rmse = torch.sqrt(mse)
    return rmse

# Testing loop
vae.eval()
total_rmse = 0
total_samples = 0

with torch.no_grad():
    for batch_idx, (data, true) in enumerate(test_loader):
        data = data.to(device)
        recon_data, _, _ = vae(data)
        rmse = calculate_rmse(recon_data, true)
        total_rmse += rmse.sum().item()
        total_samples += 32

average_rmse = total_rmse / total_samples
print(f"Average RMSE on test dataset: {average_rmse:.4f}")
