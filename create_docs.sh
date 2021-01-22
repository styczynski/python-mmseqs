#!/bin/bash

rm -rfd "./documentation" > /dev/null 2> /dev/null
rm -rfd "./__tmpdoc" > /dev/null 2> /dev/null

for dir in src/*/     # list directories in the form "/tmp/dirname/"
do
    dir=${dir%*/}      # remove the trailing "/"
    package=${dir##*/} # everything after the final "/"
    package="src/${package}"
    rm -rfd "$package/build/docs" > /dev/null 2> /dev/null
    rm -rfd "$package/pydoc-markdown.yml" > /dev/null 2> /dev/null
    poetry run bash -c "cd ${package} && pydoc-markdown --bootstrap hugo"
    poetry run bash -c "cd ${package} && pydoc-markdown"
done

cp -r .doc "./documentation"

first_src=1
for dir in src/*/
do
    dir=${dir%*/}      # remove the trailing "/"
    dir="src/${dir}"
    package=${dir##*/} # everything after the final "/"
    package_raw_name=${package}
    package="src/${package}"
    if [[ ! -f "$package/config.toml" ]]
    then
        if [[ -d "$package/build/docs" ]]
        then
            if [[ "$first_src" == "1" ]]
            then
                first_src=0
                mv "$package/build/docs/content/docs" "./documentation/content/docs/covid_genomics/$package_raw_name"
            else
                mv "$package/build/docs/content/docs" "./documentation/content/docs/covid_genomics/$package_raw_name"
            fi

            sed -i -e 1,3d "./documentation/content/docs/covid_genomics/$package_raw_name/_index.md"
            sed -i -e 1,3d "./documentation/content/docs/covid_genomics/$package_raw_name/api-documentation.md"

            printf -- '---\nweight: 10\n---' | cat - "./documentation/content/docs/covid_genomics/$package_raw_name/_index.md" > temp && mv temp "./documentation/content/docs/covid_genomics/$package_raw_name/_index.md"
            printf -- "---\ntitle: $package_raw_name\nweight: 1\n---" | cat - "./documentation/content/docs/covid_genomics/$package_raw_name/api-documentation.md" > temp && mv temp "./documentation/content/docs/covid_genomics/$package_raw_name/api-documentation.md"

            rm -f "./documentation/content/docs/covid_genomics/$package_raw_name/_index.md-e" > /dev/null 2> /dev/null
            rm -f "./documentation/content/docs/covid_genomics/$package_raw_name/api-documentation.md-e" > /dev/null 2> /dev/null

            printf -- $"- [$package_raw_name]({{< relref \"/docs/covid_genomics/$package_raw_name\" >}})\n" >> "./documentation/content/menu/index.md"
            printf -- $"  - [Usage]({{< relref \"/docs/covid_genomics/$package_raw_name/api-documentation\" >}})\n" >> "./documentation/content/menu/index.md"

        fi
    fi
done

printf -- $"<br />\n" >> "./documentation/content/menu/index.md"

#poetry run bash -c "cd amplification_model && pydoc-markdown --bootstrap hugo"