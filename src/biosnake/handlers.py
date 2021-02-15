import pandas as pd

def _client_get_databases_dataframe(self):
    data = dict(name=[], description=[], type=[])
    for db in self:
        for key in data.keys():
            data[key].append(getattr(db, key))
    return pd.DataFrame(data)

def _search_results_get_dataframe(self):
    data = dict()
    for header in self._headers:
        header_data = []
        for row in self._records:
            for aln in row:
                header_data.append(getattr(aln, header))
        data[header] = header_data
    return pd.DataFrame(data)

def _search_results_to_fasta(self, output_file_path, sequence_pattern, *header_patterns, remove_blanks=True):
    df = self.dataframe
    with open(output_file_path, 'w') as output_file:
        for index, row in df.iterrows():
            header = []
            content = row[sequence_pattern]
            if remove_blanks:
                content = content.replace('-', '')
            for header_pattern in header_patterns:
                header.append(row[header_pattern])
            output_file.write(f'>{"|".join(header)}\n')
            output_file.write(f'{content}\n')

def _database_get_dataframe(self):
    data = self.columns_data
    return pd.DataFrame(dict(
        id=data[0],
        length=data[2],
        sequence=data[1],
    ))
