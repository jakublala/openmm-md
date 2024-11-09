from notion_client import Client
import pandas as pd
import os

NOTION_TOKEN = os.getenv("NOTION_TOKEN")

# Initialize the client
notion = Client(auth=NOTION_TOKEN)

# Example: Query a database
database_id = "13922c7412678125aa83d90bc09c98b2"
response = notion.databases.query(
    database_id=database_id
)


# TODO: ideally do this, the hack is to get it from the .out file
def get_time_from_checkpoint(simulation_path):
    pass
# print(get_time_from_checkpoint("/Users/jakublala/Coding/imperial-phd/openmm-md/data/241010_FoldingUponBinding/output/241029/A-synuclein/alpha_1"))


def get_time_from_outfile(outfile_path):
    df = pd.read_csv(outfile_path, sep="\t")
    assert "(ps)" in df.columns[2], "Timestep is not in ps"
    return df["Time (ps)"].iloc[-1] / 1e-3

def populate_simulation_time():
    for result in response["results"]:
        target = result["properties"]["Target"]['select']['name'].lower()
        binder = result["properties"]["Binder"]['rich_text'][0]['plain_text'].lower()
        
        
        print(target, binder)
        # assert 0
        # simulation_path = result["properties"]["Simulation Path"]["rich_text"][0]["plain_text"]
        # time = get_time_from_outfile(os.path.join(simulation_path, "A-synuclein_alpha.out"))
        # print(time)
        # assert 0 == 1 


populate_simulation_time()









# print(len(response["results"]))

# # Example: Create a page
# notion.pages.create(
#     parent={"database_id": database_id},
#     properties={
#         "Name": {
#             "title": [
#                 {
#                     "text": {
#                         "content": "New page"
#                     }
#                 }
#             ]
#         }
#     }
# )