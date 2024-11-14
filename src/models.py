import pydantic

class Residue(pydantic.BaseModel):
    index: int
    chain_id: str
    indexing: int

class Segment(pydantic.BaseModel):
    residues: list[Residue] = pydantic.Field(default_factory=list)

    @property
    def ids(self):
        return [residue.index for residue in self.residues]

class Contact(pydantic.BaseModel):
    residue1: Residue
    residue2: Residue  

    def __str__(self):
        return f"@CA-{self.residue1.chain_id}_{self.residue1.index},@CA-{self.residue2.chain_id}_{self.residue2.index}"

class ContactMap(pydantic.BaseModel):
    contacts: list[Contact] = pydantic.Field(default_factory=list)
    
    def __str__(self):
        return '\n'.join(f"\tATOMS{i+1}={contact}" for i, contact in enumerate(self.contacts))
