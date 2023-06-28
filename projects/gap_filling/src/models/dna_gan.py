import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, Dataset
import torch.autograd as autograd
import numpy as np
from sklearn.metrics import mean_squared_error
import pandas as pd


class FruitFlyDNADataset(Dataset):
    def __init__(self, data):
        self.data = [torch.tensor(sequence, dtype=torch.float32) for sequence in data]

    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx):
        return self.data[idx]

# Generator 
class Generator(nn.Module):
    def __init__(self, input_size, hidden_size, output_size):
        super(Generator, self).__init__()
        self.net = nn.Sequential(
            nn.Linear(input_size, hidden_size*2),
            nn.LeakyReLU(0.2),
            nn.Linear(hidden_size*2, hidden_size),
            nn.LeakyReLU(0.2),
            nn.Linear(hidden_size, hidden_size),
            nn.LeakyReLU(0.2),
            nn.Linear(hidden_size, output_size),
            nn.Sigmoid()
        )

    def forward(self, x):
        return self.net(x)

# Discriminator 
class Discriminator(nn.Module):
    def __init__(self, input_size, hidden_size, output_size):
        super(Discriminator, self).__init__()
        self.net = nn.Sequential(
            nn.Linear(input_size, hidden_size),
            nn.LeakyReLU(0.2),
            nn.Linear(hidden_size, hidden_size),
            nn.LeakyReLU(0.2),
            nn.Linear(hidden_size, output_size),
            #nn.Sigmoid()
        )

    def forward(self, x):
        return self.net(x)

# Gradient penalty calculation
def compute_gradient_penalty(discriminator, real_samples, fake_samples):
    alpha = torch.Tensor(np.random.random((real_samples.size(0), 1))).to(device)
    interpolates = (alpha * real_samples + ((1 - alpha) * fake_samples)).requires_grad_(True)
    d_interpolates = discriminator(interpolates)
    gradients = autograd.grad(outputs=d_interpolates, inputs=interpolates,
                              grad_outputs=torch.ones(d_interpolates.size()).to(device),
                              create_graph=True, retain_graph=True, only_inputs=True)[0]
    gradients = gradients.view(gradients.size(0), -1)
    gradient_penalty = ((gradients.norm(2, dim=1) - 1) ** 2).mean()
    return gradient_penalty
    
# Hyperparameters
lr_g = 0.0001
lr_d = 0.001
beta1 = 0.5
beta2 = 0.9
n_critic = 1
lambda_gp = 10
input_size = 50024
hidden_size = 256
output_size = 1
batch_size = 32
epochs = 500
lr = 0.0001
k = 1

device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
print(device)

loss_function = nn.BCELoss()

data = pd.read_csv('North_America_50kMutations.tsv', delimiter='\t').iloc[:,2:].T.values
missing_data = pd.read_csv('North_America_50kMutations_5perc_missing.csv').iloc[:,3:].T.fillna(-1).values

train_data = FruitFlyDNADataset(data)
train_loader = DataLoader(train_data, batch_size=batch_size, shuffle=True)

'''
# Training loop
for epoch in range(epochs):
    for i, real_data in enumerate(train_loader):
        real_data = real_data.to(device)
        batch_size = real_data.size(0)

        # Train discriminator
        optimizer_D.zero_grad()

        real_labels = torch.ones(batch_size, 1).to(device)
        real_outputs = discriminator(real_data)
        real_loss = loss_function(real_outputs, real_labels)

        noise = torch.randn(batch_size, input_size).to(device)
        fake_data = generator(noise)
        fake_labels = torch.zeros(batch_size, 1).to(device)
        fake_outputs = discriminator(fake_data.detach())
        fake_loss = loss_function(fake_outputs, fake_labels)

        d_loss = real_loss + fake_loss
        d_loss.backward()
        optimizer_D.step()

        # Train generator
        optimizer_G.zero_grad()

        fake_outputs = discriminator(fake_data)
        g_loss = loss_function(fake_outputs, real_labels)
        g_loss.backward()
        optimizer_G.step()

    print(f"Epoch [{epoch+1}/{epochs}] - Discriminator Loss: {d_loss.item():.4f}, Generator Loss: {g_loss.item():.4f}")
'''
'''
# Training loop k:1
# 'k' number of steps to train the discriminator for each step of the generator
for epoch in range(epochs):
    for i, real_data in enumerate(train_loader):
        real_data = real_data.to(device)
        batch_size = real_data.size(0)

        # Train discriminator for k steps
        for _ in range(k):
            optimizer_D.zero_grad()

            real_labels = torch.ones(batch_size, 1).to(device)
            real_outputs = discriminator(real_data)
            real_loss = loss_function(real_outputs, real_labels)

            noise = torch.rand(batch_size, input_size).to(device)
            fake_data = generator(noise)
            fake_labels = torch.zeros(batch_size, 1).to(device)
            fake_outputs = discriminator(fake_data.detach())
            fake_loss = loss_function(fake_outputs, fake_labels)

            d_loss = real_loss + fake_loss
            d_loss.backward()
            optimizer_D.step()

        # Train generator
        optimizer_G.zero_grad()

        fake_outputs = discriminator(fake_data)
        g_loss = loss_function(fake_outputs, real_labels)
        g_loss.backward()
        optimizer_G.step()

    print(f"Epoch [{epoch+1}/{epochs}] - Discriminator Loss: {d_loss.item():.4f}, Generator Loss: {g_loss.item():.4f}")
'''
'''
# Training loop
# 'k' number of steps to train the generator for each step of the discriminator
for epoch in range(epochs):
    for i, real_data in enumerate(train_loader):
        real_data = real_data.to(device)
        batch_size = real_data.size(0)

        # Train discriminator for 1 step
        optimizer_D.zero_grad()

        real_labels = torch.ones(batch_size, 1).to(device)
        real_outputs = discriminator(real_data)
        real_loss = loss_function(real_outputs, real_labels)

        noise = torch.randn(batch_size, input_size).to(device)
        fake_data = generator(noise)
        fake_labels = torch.zeros(batch_size, 1).to(device)
        fake_outputs = discriminator(fake_data.detach())
        fake_loss = loss_function(fake_outputs, fake_labels)

        d_loss = real_loss + fake_loss
        d_loss.backward()
        optimizer_D.step()

        # Train generator for k steps
        for _ in range(k):
            optimizer_G.zero_grad()

            noise = torch.randn(batch_size, input_size).to(device)
            fake_data = generator(noise)
            fake_outputs = discriminator(fake_data)
            g_loss = loss_function(fake_outputs, real_labels)
            g_loss.backward()
            optimizer_G.step()

    print(f"Epoch [{epoch+1}/{epochs}] - Discriminator Loss: {d_loss.item():.4f}, Generator Loss: {g_loss.item():.4f}")
'''
# Initialize the models
generator = Generator(input_size, hidden_size, input_size).to(device)
discriminator = Discriminator(input_size, hidden_size, output_size).to(device)

# Define the optimizers
optimizer_G = optim.Adam(generator.parameters(), lr=lr, betas=(beta1, beta2))
optimizer_D = optim.Adam(discriminator.parameters(), lr=lr, betas=(beta1, beta2))

# Training loop
for epoch in range(epochs):
    for i, real_data in enumerate(train_loader):
        real_data = real_data.to(device)
        batch_size = real_data.size(0)

        # Train discriminator
        optimizer_D.zero_grad()

        noise = torch.randn(batch_size, input_size).to(device)
        fake_data = generator(noise)

        real_outputs = discriminator(real_data)
        fake_outputs = discriminator(fake_data.detach())

        gradient_penalty = compute_gradient_penalty(discriminator, real_data, fake_data)
        d_loss = -torch.mean(real_outputs) + torch.mean(fake_outputs) + lambda_gp * gradient_penalty
        d_loss.backward(retain_graph=True)
        optimizer_D.step()

        # Train generator every n_critic steps
        if i % n_critic == 0:
            for _ in range(k):
                optimizer_G.zero_grad()
                fake_outputs = discriminator(fake_data)
                g_loss = -torch.mean(fake_outputs)
                g_loss.backward(retain_graph=True)
                optimizer_G.step()

    print(f"Epoch [{epoch+1}/{epochs}] - Discriminator Loss: {d_loss.item():.4f}, Generator Loss: {g_loss.item():.4f}")

def impute_missing_values(sequence, generator, missing_value_indicator=-1):
    sequence = sequence.astype(np.float32)
    missing_indices = np.where(sequence == missing_value_indicator)[0]

    if len(missing_indices) > 0:
        noise = torch.rand(1, input_size).to(device)
        imputed_values = generator(noise).detach().cpu().numpy()

        for idx in missing_indices:
            sequence[idx] = imputed_values[0, idx]
    return sequence

# Impute missing values in the dataset
imputed_data = []
for i in range(missing_data.shape[0]):
    imputed_sequence = impute_missing_values(missing_data[i], generator)
    imputed_data.append(imputed_sequence)

# Calculate RMSE between imputed data and true DNA sequences
rmse = np.sqrt(mean_squared_error(imputed_data, data))
print(f"Root Mean Squared Error (RMSE): {rmse:.4f}")
