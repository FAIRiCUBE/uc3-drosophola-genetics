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

class Autoencoder(nn.Module):
    def __init__(self, num_features, hidden_dim, latent_dim):
        super(Autoencoder, self).__init__()
        self.encoder = nn.Sequential(
            nn.Linear(num_features, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, latent_dim),
            nn.ReLU()
        )
        self.decoder = nn.Sequential(
            nn.Linear(latent_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, num_features),
            nn.Sigmoid()
        )
    def forward(self, x):
        x = self.encoder(x)
        x = self.decoder(x)
        return x
        
#Hyperparameters
num_features = X_train.shape[1]
hidden_dim = 128
latent_dim = 256
epochs = 10000
model = Autoencoder(num_features, hidden_dim, latent_dim).to(device)
criterion = nn.MSELoss()
#criterion = nn.BCELoss(reduction='sum')
optimizer = optim.AdamW(model.parameters(), lr=0.001)
scheduler = MultiStepLR(optimizer, milestones=[5000], gamma=0.1)

for epoch in range(epochs):
    train_loss = 0
    model.train()
    for batch_idx, (data,) in enumerate(train_loader):
        data = data.to(device)
        optimizer.zero_grad()
        output = model(data)
        loss = criterion(output, data)
        loss.backward()
        train_loss += loss.item()
        optimizer.step()
    scheduler.step()
    print(f"Epoch [{epoch+1}/{epochs}] - Loss: {math.sqrt(train_loss / len(train_loader.dataset)):.7f}")
 
# Impute missing values
model.eval()
with torch.no_grad():
    X_test_imputed = model(X_test)
    
rmse = np.sqrt(mean_squared_error(X_train.cpu(), X_test_imputed.cpu()))
print(f'RMSE on full data: {rmse:.6f}')

X_test[mask] = X_test_imputed[mask]

missing_rmse = np.sqrt(mean_squared_error(X_train.cpu(), X_test.cpu()))
print(f'RMSE on missing data: {missing_rmse:.6f}')

X_test = np.transpose(X_test.cpu(), (1,0))

df.iloc[:, 3:] = X_test.numpy()
df.to_csv('North_America_50kMutations_5perc_missing_aefilled.csv', index=False)