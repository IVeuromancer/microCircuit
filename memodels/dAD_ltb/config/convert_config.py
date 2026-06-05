import json
import collections
def features(old_file, new_file):
    name_map = {"-50_IDRest_150":"Step_150",
                "-50_IDRest_200":"Step_200",
                "-50_IDRest_250":"Step_250",
                "-70_IDRest_150":"Step_150_hyp",
                "-50_IV_-140":"IV_-140",
                "hold_dep":"hold_dep",
                "hold_hyp":"hold_hyp",
                "RMP":"RMP"}
    old_dict = json.load(open(old_file))
    new_dict = collections.defaultdict(lambda: collections.defaultdict(list))
    for prot_name in old_dict.keys():
        for feature, values in old_dict[prot_name]["soma"].items():
            new_dict[name_map[prot_name]]["soma.v"].append({"feature":feature,
                                                            "val":values})

    with open(new_file,'w') as nw:
        json.dump(new_dict,nw, indent=1)


features("/gpfs/bbp.cscs.ch/project/proj55/iavarone/thalamus/TC_work/TC_full/config/features.json",
         "/gpfs/bbp.cscs.ch/project/proj55/software/singlecell-optimization/thalamus/config/features/bAC_TC_VPL_legacy.json")
        
