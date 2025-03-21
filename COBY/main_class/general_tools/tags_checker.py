class tags_checker:

    def tags_checker(self, dict_to_check):
        tags_list = []
        for key in ["tags", "tag"]:
            if key in dict_to_check.keys():
                tags = dict_to_check[key]
                assert type(tags) in [str, tuple, list], "Tags given to molecules must be given as either a string or a tuple or list of strings."
                if type(tags) == str:
                    tags = [tags]
                elif type(tags) == tuple:
                    tags = list(tags)
                for tag in tags:
                    assert type(tag) == str, "Tags given to molecules must be given as either a string or a tuple or list of strings."
                    assert tag not in ["tag", "tags", "all", "any"], "Tags may not be any the following words: 'tag', 'tags', 'any', 'all'."
                tags_list.extend(tags)
        return tags_list
