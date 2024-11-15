import pydantic
import numpy as np
class Residue(pydantic.BaseModel):
    index: int
    chain_id: str
    indexing: int

class Segment(pydantic.BaseModel):
    residues: list[Residue] = pydantic.Field(default_factory=list)

    @property
    def ids(self):
        return np.array([residue.index for residue in self.residues])
    
    @property
    def indexing(self):
        return self.residues[0].indexing

    # validation
    @pydantic.validate_arguments
    def validate(self):
        assert np.all(residue.indexing == self.indexing for residue in self.residues), "All residues must have the same indexing"

    def __str__(self):
        return ','.join(f"@CA-{residue.chain_id}_{residue.index}" for residue in self.residues)

class Contact(pydantic.BaseModel):
    residue1: Residue
    residue2: Residue  

    def __str__(self):
        return f"@CA-{self.residue1.chain_id}_{self.residue1.index},@CA-{self.residue2.chain_id}_{self.residue2.index}"

class ContactMap(pydantic.BaseModel):
    contacts: list[Contact] = pydantic.Field(default_factory=list)
    
    def __str__(self):
        return '\n'.join(f"\tATOMS{i+1}={contact}" for i, contact in enumerate(self.contacts))
