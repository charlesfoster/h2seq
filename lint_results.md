## `nf-core pipelines lint` overall result: Failed :x:

Posted for pipeline commit 2794dc9

```diff
+| ✅ 192 tests passed       |+
#| ❔  28 tests were ignored |#
!| ❗   9 tests had warnings |!
-| ❌   1 tests failed       |-
```

<details>

### :x: Test failures:

- [multiqc_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/multiqc_config) - `assets/multiqc_config.yml` does not meet requirements: Section charlesfoster-h2seq-summary should have the lowest order

### :heavy_exclamation_mark: Test warnings:

- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config `manifest.version` should end in `dev`: `0.1.2`
- [readme](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/readme) - README did not have an nf-core template version badge.
- [pipeline_todos](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/pipeline_todos) - TODO string in `main.nf`: _add in other checks here based on all params_
- [pipeline_todos](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/pipeline_todos) - TODO string in `main.nf`: _Optionally add in-text citation tools to this list._
- [pipeline_todos](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/pipeline_todos) - TODO string in `main.nf`: _Optionally add bibliographic entries to this list._
- [pipeline_todos](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/pipeline_todos) - TODO string in `main.nf`: _Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!_
- [local_component_structure](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/local_component_structure) - long_read_mapping.nf in subworkflows/local should be moved to a SUBWORKFLOW_NAME/main.nf structure
- [local_component_structure](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/local_component_structure) - prepare_consensus_bam.nf in subworkflows/local should be moved to a SUBWORKFLOW_NAME/main.nf structure
- [local_component_structure](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/local_component_structure) - short_read_mapping.nf in subworkflows/local should be moved to a SUBWORKFLOW_NAME/main.nf structure

### :grey_question: Tests ignored:

- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File is ignored: `CODE_OF_CONDUCT.md`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File is ignored: `assets/nf-core-h2seq_logo_light.png`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File is ignored: `docs/images/nf-core-h2seq_logo_light.png`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File is ignored: `docs/images/nf-core-h2seq_logo_dark.png`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File is ignored: `.github/ISSUE_TEMPLATE/config.yml`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File is ignored: `.github/workflows/awstest.yml`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File is ignored: `.github/workflows/awsfulltest.yml`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File is ignored: `conf/igenomes.config`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File is ignored: `.github/workflows/nf-test.yml`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File is ignored: `.github/actions/get-shards/action.yml`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File is ignored: `.github/actions/nf-test/action.yml`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File is ignored: `nf-test.config`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File is ignored: `tests/default.nf.test`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File is ignored: `.github/workflows/release-announcements.yml`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File is ignored: `conf/igenomes_ignored.config`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File is ignored: `ro-crate-metadata.json`
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config variable ignored: `manifest.name`
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config variable ignored: `manifest.version`
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config variable ignored: `manifest.homePage`
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config variable ignored: `process.cpus`
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config variable ignored: `process.memory`
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config variable ignored: `process.time`
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config variable ignored: `params.max_cpus`
- [nf_test_content](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nf_test_content) - nf_test_content
- [files_unchanged](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_unchanged) - Required pipeline config not found - {'manifest.contributors'}
- [actions_nf_test](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/actions_nf_test) - '.github/workflows/nf-test.yml' not found
- [actions_awstest](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/actions_awstest) - 'awstest.yml' workflow not found: `.github/workflows/awstest.yml`
- [rocrate_readme_sync](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/rocrate_readme_sync) - `ro-crate-metadata.json` not found

### :white_check_mark: Tests passed:

- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `.gitattributes`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `.gitignore`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `.nf-core.yml`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `.prettierignore`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `.prettierrc.yml`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `CHANGELOG.md`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `CITATIONS.md`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `LICENSE` or `LICENSE.md` or `LICENCE` or `LICENCE.md`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `nextflow_schema.json`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `nextflow.config`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `README.md`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `.github/.dockstore.yml`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `.github/CONTRIBUTING.md`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `.github/ISSUE_TEMPLATE/bug_report.yml`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `.github/ISSUE_TEMPLATE/feature_request.yml`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `.github/PULL_REQUEST_TEMPLATE.md`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `.github/workflows/branch.yml`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `.github/workflows/linting_comment.yml`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `.github/workflows/linting.yml`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `assets/email_template.html`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `assets/email_template.txt`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `assets/sendmail_template.txt`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `conf/modules.config`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `conf/test.config`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `conf/test_full.config`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `docs/output.md`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `docs/README.md`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `docs/README.md`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `docs/usage.md`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `main.nf`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `assets/multiqc_config.yml`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `conf/base.config`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File found: `modules.json`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File not found check: `.github/ISSUE_TEMPLATE/bug_report.md`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File not found check: `.github/ISSUE_TEMPLATE/feature_request.md`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File not found check: `.github/workflows/push_dockerhub.yml`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File not found check: `.markdownlint.yml`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File not found check: `.nf-core.yaml`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File not found check: `.yamllint.yml`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File not found check: `bin/markdown_to_html.r`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File not found check: `conf/aws.config`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File not found check: `docs/images/nf-core-h2seq_logo.png`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File not found check: `lib/Checks.groovy`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File not found check: `lib/Completion.groovy`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File not found check: `lib/NfcoreTemplate.groovy`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File not found check: `lib/Utils.groovy`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File not found check: `lib/Workflow.groovy`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File not found check: `lib/WorkflowMain.groovy`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File not found check: `lib/WorkflowH2seq.groovy`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File not found check: `parameters.settings.json`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File not found check: `pipeline_template.yml`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File not found check: `Singularity`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File not found check: `lib/nfcore_external_java_deps.jar`
- [files_exist](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/files_exist) - File not found check: `.travis.yml`
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Found nf-schema plugin
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config variable found: `manifest.nextflowVersion`
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config variable found: `manifest.description`
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config variable found: `timeline.enabled`
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config variable found: `trace.enabled`
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config variable found: `report.enabled`
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config variable found: `dag.enabled`
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config variable found: `params.outdir`
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config variable found: `params.input`
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config variable found: `manifest.mainScript`
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config variable found: `timeline.file`
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config variable found: `trace.file`
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config variable found: `report.file`
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config variable found: `dag.file`
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config variable (correctly) not found: `params.nf_required_version`
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config variable (correctly) not found: `params.container`
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config variable (correctly) not found: `params.singleEnd`
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config variable (correctly) not found: `params.igenomesIgnore`
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config variable (correctly) not found: `params.name`
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config variable (correctly) not found: `params.enable_conda`
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config `timeline.enabled` had correct value: `true`
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config `report.enabled` had correct value: `true`
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config `trace.enabled` had correct value: `true`
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config `dag.enabled` had correct value: `true`
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config `dag.file` ended with `.html`
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config variable `manifest.nextflowVersion` started with >= or !>=
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - nextflow.config contains configuration profile `test`
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config default value correct: params.outdir= results
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config default value correct: params.long_reads_min_len= 500
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config default value correct: params.long_reads_max_len= 2000
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config default value correct: params.long_reads_min_qual= 9.0
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config default value correct: params.short_reads_min_len= 70
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config default value correct: params.fastp_qualified_quality= 15
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config default value correct: params.fastp_unqualified_percent_limit= 40
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config default value correct: params.fastp_cut_mean_quality= 15
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config default value correct: params.fastp_low_complexity_filter= true
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config default value correct: params.fastp_complexity_threshold= 40
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config default value correct: params.fastp_error_correction= true
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config default value correct: params.reference_selection_tool= salmon
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config default value correct: params.abundance_top_percentage= 5.0
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config default value correct: params.clip_tolerance= 5
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config default value correct: params.consensus_min_depth= 15
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config default value correct: params.ont_min_snv_af= 0.15
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config default value correct: params.ont_min_indel_af= 0.5
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config default value correct: params.illumina_min_snv_af= 0.05
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config default value correct: params.illumina_min_indel_af= 0.5
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config default value correct: params.consensus_call_af= 0.5
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config default value correct: params.majority_allele_consensus= false
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config default value correct: params.ont_variant_min_qual= 5.0
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config default value correct: params.illumina_variant_min_qual= 20.0
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config default value correct: params.clair3_model_path= /opt/models/r1041_e82_400bps_sup_v410
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config default value correct: params.clair3_platform= ont
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config default value correct: params.clair3_chunk_size= 15000
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config default value correct: params.custom_config_version= master
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config default value correct: params.custom_config_base= https://raw.githubusercontent.com/nf-core/configs/master
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config default value correct: params.max_cpus= 16
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config default value correct: params.max_memory= 60.GB
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config default value correct: params.max_time= 240.h
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config default value correct: params.publish_dir_mode= copy
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config default value correct: params.max_multiqc_email_size= 25.MB
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config default value correct: params.validate_params= true
- [nextflow_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nextflow_config) - Config default value correct: params.pipelines_testdata_base_path= https://github.com/charlesfoster/h2seq/tree/dev/test-datasets
- [readme](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/readme) - README Nextflow minimum version badge matched config. Badge: `23.10.0`, Config: `23.10.0`
- [readme](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/readme) - README Zenodo placeholder was replaced with DOI.
- [pipeline_if_empty_null](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/pipeline_if_empty_null) - No `ifEmpty(null)` strings found
- [plugin_includes](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/plugin_includes) - No wrong validation plugin imports have been found
- [pipeline_name_conventions](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/pipeline_name_conventions) - Name adheres to nf-core convention
- [template_strings](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/template_strings) - Did not find any Jinja template strings (0 files)
- [schema_lint](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/schema_lint) - Schema lint passed
- [schema_lint](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/schema_lint) - Schema title + description lint passed
- [schema_lint](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/schema_lint) - Input mimetype lint passed: 'text/csv'
- [schema_params](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/schema_params) - Schema matched params returned from nextflow config
- [system_exit](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/system_exit) - No `System.exit` calls found
- [actions_schema_validation](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/actions_schema_validation) - Workflow validation passed: linting_comment.yml
- [actions_schema_validation](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/actions_schema_validation) - Workflow validation passed: branch.yml
- [actions_schema_validation](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/actions_schema_validation) - Workflow validation passed: fix-linting.yml
- [actions_schema_validation](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/actions_schema_validation) - Workflow validation passed: clean-up.yml
- [actions_schema_validation](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/actions_schema_validation) - Workflow validation passed: linting.yml
- [actions_schema_validation](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/actions_schema_validation) - Workflow validation passed: ci.yml
- [actions_schema_validation](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/actions_schema_validation) - Workflow validation passed: download_pipeline.yml
- [merge_markers](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/merge_markers) - No merge markers found in pipeline files
- [modules_json](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/modules_json) - Only installed modules found in `modules.json`
- [multiqc_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/multiqc_config) - `assets/multiqc_config.yml` found and not ignored.
- [multiqc_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/multiqc_config) - `assets/multiqc_config.yml` contains `report_section_order`
- [multiqc_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/multiqc_config) - `assets/multiqc_config.yml` contains `export_plots`
- [multiqc_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/multiqc_config) - `assets/multiqc_config.yml` contains `report_comment`
- [multiqc_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/multiqc_config) - `assets/multiqc_config.yml` contains 'export_plots: true'.
- [modules_structure](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/modules_structure) - modules directory structure is correct 'modules/nf-core/TOOL/SUBTOOL'
- [local_component_structure](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/local_component_structure) - local modules directory structure is correct 'modules/local/TOOL/SUBTOOL'
- [base_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/base_config) - `conf/base.config` found and not ignored.
- [modules_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/modules_config) - `conf/modules.config` found and not ignored.
- [modules_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/modules_config) - `FASTQC` found in `conf/modules.config` and Nextflow scripts.
- [modules_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/modules_config) - `MULTIQC` found in `conf/modules.config` and Nextflow scripts.
- [modules_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/modules_config) - `BUILD_MULTIQC_SECTIONS` found in `conf/modules.config` and Nextflow scripts.
- [modules_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/modules_config) - `NANOQ` found in `conf/modules.config` and Nextflow scripts.
- [modules_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/modules_config) - `FASTP` found in `conf/modules.config` and Nextflow scripts.
- [modules_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/modules_config) - `SELECT_BEST_REFERENCE` found in `conf/modules.config` and Nextflow scripts.
- [modules_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/modules_config) - `REFERENCE_METADATA_FROM_FASTA` found in `conf/modules.config` and Nextflow scripts.
- [modules_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/modules_config) - `SEQKIT_GREP` found in `conf/modules.config` and Nextflow scripts.
- [modules_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/modules_config) - `CALCULATE_READ_STATS` found in `conf/modules.config` and Nextflow scripts.
- [modules_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/modules_config) - `SEQKIT_STATS` found in `conf/modules.config` and Nextflow scripts.
- [modules_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/modules_config) - `SEQKIT_STATS_RAW_LONG` found in `conf/modules.config` and Nextflow scripts.
- [modules_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/modules_config) - `SEQKIT_STATS_CLEAN_LONG` found in `conf/modules.config` and Nextflow scripts.
- [modules_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/modules_config) - `SEQKIT_STATS_RAW_SHORT` found in `conf/modules.config` and Nextflow scripts.
- [modules_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/modules_config) - `SEQKIT_STATS_CLEAN_SHORT` found in `conf/modules.config` and Nextflow scripts.
- [modules_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/modules_config) - `KALLISTO_INDEX` found in `conf/modules.config` and Nextflow scripts.
- [modules_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/modules_config) - `KALLISTO_QUANT` found in `conf/modules.config` and Nextflow scripts.
- [modules_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/modules_config) - `SALMON_INDEX` found in `conf/modules.config` and Nextflow scripts.
- [modules_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/modules_config) - `SALMON_QUANT` found in `conf/modules.config` and Nextflow scripts.
- [modules_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/modules_config) - `BWA_INDEX` found in `conf/modules.config` and Nextflow scripts.
- [modules_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/modules_config) - `BWA_MEM` found in `conf/modules.config` and Nextflow scripts.
- [modules_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/modules_config) - `BEDTOOLS_BAMTOBED` found in `conf/modules.config` and Nextflow scripts.
- [modules_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/modules_config) - `SPLIT_CONSENSUS_GENOMES` found in `conf/modules.config` and Nextflow scripts.
- [modules_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/modules_config) - `SAMTOOLS_AMPLICONCLIP` found in `conf/modules.config` and Nextflow scripts.
- [modules_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/modules_config) - `MINIMAP2_ALIGN` found in `conf/modules.config` and Nextflow scripts.
- [modules_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/modules_config) - `SAMTOOLS_SORT` found in `conf/modules.config` and Nextflow scripts.
- [modules_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/modules_config) - `CLAIR3` found in `conf/modules.config` and Nextflow scripts.
- [modules_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/modules_config) - `LOFREQ_INDELQUAL` found in `conf/modules.config` and Nextflow scripts.
- [modules_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/modules_config) - `LOFREQ_CALL` found in `conf/modules.config` and Nextflow scripts.
- [modules_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/modules_config) - `PREPARE_CLAIR3_VCF` found in `conf/modules.config` and Nextflow scripts.
- [modules_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/modules_config) - `PREPARE_LOFREQ_VCF` found in `conf/modules.config` and Nextflow scripts.
- [modules_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/modules_config) - `FILTER_VARIANTS` found in `conf/modules.config` and Nextflow scripts.
- [modules_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/modules_config) - `ANNOTATE_VARIANTS` found in `conf/modules.config` and Nextflow scripts.
- [modules_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/modules_config) - `CREATE_CONSENSUS_MASK` found in `conf/modules.config` and Nextflow scripts.
- [modules_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/modules_config) - `BCFTOOLS_CONSENSUS` found in `conf/modules.config` and Nextflow scripts.
- [modules_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/modules_config) - `REMOVE_EMPTY_SEQUENCES` found in `conf/modules.config` and Nextflow scripts.
- [modules_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/modules_config) - `HCV_GLUE` found in `conf/modules.config` and Nextflow scripts.
- [modules_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/modules_config) - `GENERATE_WHOLE_GENOME_BED` found in `conf/modules.config` and Nextflow scripts.
- [modules_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/modules_config) - `MOSDEPTH_GENOME` found in `conf/modules.config` and Nextflow scripts.
- [modules_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/modules_config) - `COVERAGE_METRICS` found in `conf/modules.config` and Nextflow scripts.
- [modules_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/modules_config) - `PARSE_HCV_GLUE_COVERAGE` found in `conf/modules.config` and Nextflow scripts.
- [modules_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/modules_config) - `PLOT_HCV_SUMMARY` found in `conf/modules.config` and Nextflow scripts.
- [modules_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/modules_config) - `PLOT_DEPTH_SUMMARY` found in `conf/modules.config` and Nextflow scripts.
- [modules_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/modules_config) - `SAMTOOLS_FAIDX` found in `conf/modules.config` and Nextflow scripts.
- [modules_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/modules_config) - `COUNT_MAPPED_READS` found in `conf/modules.config` and Nextflow scripts.
- [modules_config](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/modules_config) - `BUILD_RUN_SUMMARY` found in `conf/modules.config` and Nextflow scripts.
- [nfcore_yml](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nfcore_yml) - Repository type in `.nf-core.yml` is valid: `pipeline`
- [nfcore_yml](https://nf-co.re/tools/docs/3.5.2/pipeline_lint_tests/nfcore_yml) - nf-core version in `.nf-core.yml` is set to the latest version: `3.5.2`

### Run details

- nf-core/tools version 3.5.2
- Run at `2026-03-18 14:07:45`

</details>
