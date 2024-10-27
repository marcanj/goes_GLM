import eumdac
import datetime
import shutil

# Define EUMETSAT API credentials
consumer_key = "EZlHthGQ2GqUVy2t5avtxOpUdMAa"  # Replace with your Consumer Key
consumer_secret = "ezIF5cgknJWQixdBuIxkfkcaSDEa"  # Replace with your Consumer Secret

# Feed the token object with your credentials, find yours at https://api.eumetsat.int/api-key/
credentials = (consumer_key, consumer_secret)
token = eumdac.AccessToken(credentials)

# Create datastore object with with your token
datastore = eumdac.DataStore(token)
# Select an FCI collection, eg "FCI Level 1c High Resolution Image Data - MTG - 0 degree" - "EO:EUM:DAT:0665"
selected_collection = datastore.get_collection('EO:EUM:DAT:0691')

# Set sensing start and end time
start = datetime.datetime(2024, 10, 5, 00, 00)
end = datetime.datetime(2024, 10, 5, 00, 5)

# Retrieve datasets that match the filter
products = selected_collection.search(
    dtstart=start,
    dtend=end)

# Print found products
for product in products:
        print(product)

# Download all found products
for product in products:
    with product.open() as source_file, open(source_file.name, mode='wb') as destination_file:
        shutil.copyfileobj(source_file, destination_file)
        print(f'Download of product {product} finished.')