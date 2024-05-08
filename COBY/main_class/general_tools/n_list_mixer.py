import itertools

class n_list_mixer:
    def n_list_mixer(self, *Lists):
        ### Mixes n lists to achieve even spacing between values from initial lists
        mix_determiner = []
        for l in Lists:
            mix_determiner.append([(i / len(l), v) for i, v in enumerate(l)])
        mixed = [v for i, v in sorted(list(itertools.chain(*mix_determiner)))]
        return mixed

