import reve

# filename = 'inputs/SiO2/15120/sio2-15120at-1frame/default_1frame.xyz'
filename = 'inputs/SiO2/15120/sio2-15120at-51_frames/default.xyz'
# filename = 'inputs/SiO2/27216/sample-27216at-100frames-npt-0GPa/sample-27216at-100frames-npt-0GPa.lammpstrj'

settings = (reve.settings.SettingsBuilder() \
    .with_project_name('test') \
    .with_export_directory('export') \
    .with_file_location(filename) \
    .with_number_of_atoms(15120) \
    .with_range_of_frames(0, -1) \
    .with_quiet(False) \
    .build()
)

reve.main(settings)
