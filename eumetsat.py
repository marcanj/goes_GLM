import requests

# Define EUMETSAT API credentials
consumer_key = "EZlHthGQ2GqUVy2t5avtxOpUdMAa"  # Replace with your Consumer Key
consumer_secret = "ezIF5cgknJWQixdBuIxkfkcaSDEa"  # Replace with your Consumer Secret

# Encode the authorization header (Basic auth)
import base64
auth_string = f"{consumer_key}:{consumer_secret}"
auth_header = base64.b64encode(auth_string.encode()).decode()

# Step 1: Obtain the API token
token_url = "https://api.eumetsat.int/token"
token_headers = {
    "Authorization": f"Basic {auth_header}"
}
token_data = {
    "grant_type": "client_credentials"
}
token_response = requests.post(token_url, headers=token_headers, data=token_data, verify=False)
if token_response.status_code == 200:
    access_token = token_response.json()["access_token"]
    print("Token retrieved successfully:", access_token)
else:
    print("Failed to retrieve token:", token_response.status_code, token_response.text)
    exit()

# Step 2: Define headers for further API requests with the token
api_headers = {
    "Authorization": f"Bearer {access_token}"
}

# Define API endpoints
data_search_url = "https://api.eumetsat.int/data-store/search"
data_order_url = "https://api.eumetsat.int/data-store/order"
data_download_url = "https://api.eumetsat.int/data-store/download"

# Step 3: Search for Lightning Imager data
search_params = {
    "query": "Lightning Imager",  # Update search parameters if needed
    "collectionId": "EO:EUM:DAT:MTG:LI",  # Example collection ID for LI data
    "startDate": "2024-10-05T00:00:00Z",
    "endDate": "2024-10-06T00:00:00Z",
    "bbox": "-10,35,5,45"  # Define bounding box for area of interest
}
search_response = requests.get(data_search_url, headers=api_headers, params=search_params)
if search_response.status_code == 200:
    datasets = search_response.json()
    print("Datasets found:", datasets)
    dataset_id = datasets[0]["id"]  # Select the first dataset ID (or refine as needed)
else:
    print("Dataset search failed:", search_response.status_code)
    exit()

# Step 4: Place an order for the dataset
order_payload = {
    "datasetId": dataset_id,
    "format": "application/x-netcdf"  # Change if needed
}
order_response = requests.post(data_order_url, headers=api_headers, json=order_payload)
if order_response.status_code == 200:
    order_id = order_response.json()["id"]
    print("Order placed successfully, ID:", order_id)
else:
    print("Order placement failed:", order_response.status_code)
    exit()

# Step 5: Check the order status and download if ready
status_url = f"{data_order_url}/{order_id}/status"
status_response = requests.get(status_url, headers=api_headers)
if status_response.status_code == 200 and status_response.json()["status"] == "completed":
    download_response = requests.get(f"{data_download_url}/{order_id}", headers=api_headers)
    filename = "lightning_data.nc"
    with open(filename, "wb") as file:
        file.write(download_response.content)
    print(f"Data downloaded and saved as {filename}")
else:
    print("Order is not ready, status:", status_response.json()["status"])
