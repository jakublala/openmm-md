from notion_client import Client
import pandas as pd
import os
from tqdm import tqdm
from dotenv import load_dotenv
load_dotenv()

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

target_mapping = {
    "CD28": "CD28",
    "A-synuclein": "ASYN",
    "p53": "P53",
    "SUMO": "SUMO"
}
target_mapping = {k.lower(): v for k, v in target_mapping.items()}

binder_mapping = {
    "alpha": "A",
    "beta": "B",
    "general": "G",
    "partial": "P",
    "end": "E",
    "1": "1",
    "2": "2",
    "1c": "1c",
    "1a": "1a"
}
binder_mapping = {k.lower(): v for k, v in binder_mapping.items()}

def get_time_from_outfile(outfile_path):
    df = pd.read_csv(outfile_path, sep="\t")
    assert "(ps)" in df.columns[2], "Timestep is not in ps"
    return round(df["Time (ps)"].iloc[-1] / 1e3, 0)

def assign_time(id, time, property_name="Time [ns]"):
    try:
        notion.pages.update(
            page_id=id,
            properties={property_name: time}
        )
        return True
    except Exception as e:
        print(f"Error assigning time to {id}: {e}")
        return False


def populate_trajectories():
    pass

def get_id(target, binder):
    return f"{target_mapping[target]}-{binder_mapping[binder]}"

def upload_trajectory(id, filepath):
    url = filepath.replace("../../", "", 1)
    # add the
    # https://imperiallondon-my.sharepoint.com/:f:/r/personal/jl24018_ic_ac_uk/Documents/data?csf=1&web=1&e=RJG5En
    print(f"Uploading {url}")
    try:
        notion.pages.upload_file(
            file=filepath,
            parent_id=id
        )
        return True
    except Exception as e:
        print(f"Error uploading trajectory to {id}: {e}")
        return False

num_runs = 5
def populate_simulation_time():
    
    for result in tqdm(response["results"], desc="Processing results"):
        target = result["properties"]["Target"]['select']['name']
        binder = result["properties"]["Binder"]['rich_text'][0]['plain_text']
        date = result["properties"]["Date"]['number']
        page_id = result["id"]
        title = result["properties"]["ID"]['title'][0]["plain_text"]

        num_run = int(title.split("-")[-1])

        if target.lower() in ["p53", "sumo"]:
            continue
        
        folder = f"../../data/241010_FoldingUponBinding/output/{date}/{target}/{binder}_{num_run}"
        # assert tthe folder exists and it has a .out file
        assert os.path.exists(folder), f"Folder {folder} does not exist"
        assert os.path.exists(os.path.join(folder, f"{target}_{binder}.out")), f"File {binder}.out does not exist in {folder}"
        time = get_time_from_outfile(os.path.join(folder, f"{target}_{binder}.out"))
        # assign_time(page_id, time)
        upload_trajectory(page_id, os.path.join(folder, f"{target}_{binder}.dcd"))

        
        
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