#!/bin/bash
set -e

# Configuration
RESOURCE_GROUP="rg-md-eastus"
LOCATION="eastus"
STORAGE_ACCOUNT_NAME="stmdsimulations$RANDOM" # Random suffix to ensure uniqueness
CONTAINER_NAME="simulations"

echo "=== Setting up Azure Storage for Simulation Results ==="

# 1. Create Storage Account
echo "Creating Storage Account: $STORAGE_ACCOUNT_NAME in $RESOURCE_GROUP..."
az storage account create \
    --name $STORAGE_ACCOUNT_NAME \
    --resource-group $RESOURCE_GROUP \
    --location $LOCATION \
    --sku Standard_LRS \
    --kind StorageV2

# 2. Create Container
echo "Creating Container: $CONTAINER_NAME..."
# We need the account key to create the container or just use auth mode login if role assigned. 
# Simplest is grabbing connection string or key.
AZURE_STORAGE_KEY=$(az storage account keys list --resource-group $RESOURCE_GROUP --account-name $STORAGE_ACCOUNT_NAME --query "[0].value" -o tsv)

az storage container create \
    --name $CONTAINER_NAME \
    --account-name $STORAGE_ACCOUNT_NAME \
    --account-key $AZURE_STORAGE_KEY

echo ""
echo "=== Storage Setup Complete ==="
echo "Storage Account: $STORAGE_ACCOUNT_NAME"
echo "Container:       $CONTAINER_NAME"
echo ""
echo "Please add the following line to your .env file:"
echo "STORAGE_ACCOUNT_NAME=$STORAGE_ACCOUNT_NAME"
echo "STORAGE_ACCOUNT_KEY=$AZURE_STORAGE_KEY"
