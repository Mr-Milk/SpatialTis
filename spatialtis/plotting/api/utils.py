def query_df(data, query):
    query_statement = []
    for k, v in query.items():
        if isinstance(v, str):
            query_statement.append(f"({k}=='{v}')")
        else:
            query_statement.append(f"({k}=={v})")
    return data.query("&".join(query_statement)).copy()
