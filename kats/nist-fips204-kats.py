import json

def write_mldsa_kat_text(json_file, output_file, tg_id):
    with open(json_file, "r") as f:
        data = json.load(f)

    with open(output_file, "w") as out:
        for group in data.get("testGroups", []):
            if group.get("tgId") == tg_id:
                for test in group.get("tests", []):
                    out.write(f"tid = {test['tcId']:02}\n")
                    out.write(f"xi = {test['seed']}\n")
                    out.write(f"pk = {test['pk']}\n")
                    out.write(f"sk = {test['sk']}\n\n")

if __name__ == "__main__":
    json_file = "acvp-ml-dsa-keygen-fips204.json"
    output_file = "nist-acvp-keygen44-kats.txt"
    write_mldsa_kat_text(json_file, output_file, 1)
    print(f"Written output to {output_file}")
    output_file = "nist-acvp-keygen65-kats.txt"
    write_mldsa_kat_text(json_file, output_file, 2)
    print(f"Written output to {output_file}")
    output_file = "nist-acvp-keygen87-kats.txt"
    write_mldsa_kat_text(json_file, output_file, 3)
    print(f"Written output to {output_file}")
