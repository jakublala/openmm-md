import pydantic

# TODO: create scheme here that will then turn themselves into a plumed.dat file
# have some defaults but then also require a lot of the values
# extract them from some initial config file.
# use hydra-core to do this
# https://hydra.cc/docs/intro
class PlumedConfigSchema(pydantic.BaseModel):
    type: str


class ContactMapSchema(pydantic.BaseModel):
    cutoff: float
    contact_residues: list[int]

