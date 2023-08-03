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

class CNNAutoencoder(nn.Module):
    def __init__(self):
        super(CNNAutoencoder, self).__init__()
        self.encoder = nn.Sequential(
            nn.Conv1d(1, 32, kernel_size=3, stride=2, padding=0), 
            nn.ReLU(),
            nn.Conv1d(32, 64, kernel_size=3, stride=2, padding=0), 
            nn.ReLU(),
            nn.Conv1d(64, 128, kernel_size=3, stride=2, padding=0), 
            nn.ReLU(),
            nn.Conv1d(128, 256, kernel_size=3, stride=2, padding=0), 
            nn.ReLU(),
            nn.Conv1d(256, 512, kernel_size=3, stride=2, padding=0), 
            nn.ReLU(),
            nn.Conv1d(512, 1024, kernel_size=3, stride=2, padding=0), 
            nn.ReLU(),
            nn.Conv1d(1024, 2048, kernel_size=3, stride=2, padding=0), 
            nn.ReLU(),
        )
        self.decoder = nn.Sequential(
            nn.ConvTranspose1d(2048, 1024, kernel_size=3, stride=2, padding=0, output_padding=1), 
            nn.ReLU(),
            nn.ConvTranspose1d(1024, 512, kernel_size=3, stride=2, padding=0, output_padding=1), 
            nn.ReLU(),
            nn.ConvTranspose1d(512, 256, kernel_size=3, stride=2, padding=0, output_padding=0), 
            nn.ReLU(),
            nn.ConvTranspose1d(256, 128, kernel_size=3, stride=2, padding=0, output_padding=1), 
            nn.ReLU(),
            nn.ConvTranspose1d(128, 64, kernel_size=3, stride=2, padding=0, output_padding=0), 
            nn.ReLU(),
            nn.ConvTranspose1d(64, 32, kernel_size=3, stride=2, padding=0, output_padding=1), 
            nn.ReLU(),
            nn.ConvTranspose1d(32, 1, kernel_size=3, stride=2, padding=0, output_padding=0),
        )
        self.sigmoid = nn.Sigmoid()

    def forward(self, x):
        x = x.unsqueeze(1)
        x = self.encoder(x)
        #print(x.size())
        x = self.decoder(x)
        #print(x.size())
        x = x.squeeze(1)
        return self.sigmoid(x)[:,:-1]

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
epochs = 1000
model = CNNAutoencoder().to(device)
#criterion = nn.MSELoss()
criterion = nn.BCELoss(reduction='sum')
#criterion = nn.MSELoss(reduction='none')
optimizer = optim.Adam(model.parameters(), lr=0.0001)
scheduler = MultiStepLR(optimizer, milestones=[500], gamma=0.1)
'''
# Train the model
for epoch in range(epochs):
    model.train()
    optimizer.zero_grad()
    output = model(X_train)
    loss = criterion(output, X_train)
    loss = loss[~mask]  
    loss = loss.mean() 
    loss.backward()
    optimizer.step()
    #scheduler.step()
    print(f'epoch {epoch+1}/{epochs}, loss={loss.item():.7f}')
'''
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
    #scheduler.step()
    print(f"Epoch [{epoch+1}/{epochs}] - Loss: {math.sqrt(train_loss / len(train_loader.dataset)):.7f}")
    
# Impute missing values
model.eval()
with torch.no_grad():
    X_test_imputed = model(X_test)

rmse = math.sqrt(mean_squared_error(X_train.cpu(), X_test_imputed.cpu()))
print(f'RMSE on full data: {rmse:.6f}')

X_test[mask] = X_test_imputed[mask]

missing_rmse = math.sqrt(mean_squared_error(X_train.cpu(), X_test.cpu()))
print(f'RMSE on missing data: {missing_rmse:.6f}')

#np.save('North_America_50kMutations_5perc_missing_maefilled.npy', X_test.cpu())

X_test = np.transpose(X_test.cpu(), (1,0))

df.iloc[:, 3:] = X_test.numpy()
df.to_csv('North_America_50kMutations_5perc_missing_cnnfilled.csv', index=False)