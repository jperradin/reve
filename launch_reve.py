import reve

settings = (reve.settings.SettingsBuilder() \
    .with_project_name('test') \
    .with_export_directory('export') \
    .with_file_location('./inputs/51_frames/default.xyz') \
    .with_number_of_atoms(15120) \
    .with_quiet(False) \
    .with_overwrite_results(False) \
    .build()
)

reve.main(settings)
