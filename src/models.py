import pydantic
import numpy as np
from typing import Optional, List

class Residue(pydantic.BaseModel):
    index: int
    global_index: int
    chain_id: str
    indexing: int
    atom_indices: Optional[List[int]] = None

    def __str__(self):
        return f"@CA-{self.chain_id}_{self.index}"

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
    
    @pydantic.validate_arguments
    def validate_length(self):
        assert len(self.contacts) > 0, "Contact map must contain at least one contact"
    
    @pydantic.validate_arguments
    def validate_chain_ids(self):
        residues1 = self.all_residues1
        residues2 = self.all_residues2
        assert all(res.chain_id == residues1[0].chain_id for res in residues1), "All residues in residues1 must have the same chain ID"
        assert all(res.chain_id == residues2[0].chain_id for res in residues2), "All residues in residues2 must have the same chain ID"
    

    @property
    def all_residues1(self):
        return [contact.residue1 for contact in self.contacts]
    
    @property
    def all_residues2(self):
        return [contact.residue2 for contact in self.contacts]