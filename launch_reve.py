import reve

settings = reve.settings.Settings()

settings.project_name.set_value('test')
settings.export_directory.set_value('export')
settings.file_locations.set_value('./inputs/51_frames/')
settings.file_format.set_value('xyz')
settings.composition.set_value('SiO2')
settings.number_of_frames.set_value(51)
settings.number_of_atoms.set_value(15120)
settings.range_of_frames.set_value((0, 50))
settings.timestep.set_value(0.0016)
settings.temperature.set_value(300)
settings.pressure.set_value(0)
settings.volume.set_value(0)
settings.quiet.set_value(False)
settings.overwrite_results.set_value(False)

analyzers = [
    'pair_distribution_function',
    # 'bond_angular_distribution',
    # 'neutron_structure_factor'
]

settings.analyzers.set_value(analyzers)

reve.main(settings)

