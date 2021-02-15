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

def resolve_vars(text, keys, vars):
    for key in keys:
        key_val = '{'+key+'}'
        if key_val in text:
            text = text.replace('{'+key+'}', str(vars[key]))
    return text

def _search_results_to_fasta(self, output_file_path, header_pattern, sequence_pattern, remove_blanks=True):
    df = self.dataframe
    outputs = dict()
    for index, row in df.iterrows():
        header = resolve_vars(header_pattern, df.columns, row)
        content = resolve_vars(sequence_pattern, df.columns, row)
        output_path = resolve_vars(output_file_path, df.columns, row)
        if output_path not in outputs:
            outputs[output_path] = open(output_path, 'w')
        if remove_blanks:
            content = content.replace('-', '')
        outputs[output_path].write(f'>{header}\n')
        outputs[output_path].write(f'{content}\n')
    for output_file in outputs.values():
        output_file.close()

def _database_get_dataframe(self):
    data = self.columns_data
    return pd.DataFrame(dict(
        id=data[0],
        length=data[2],
        sequence=data[1],
    ))
