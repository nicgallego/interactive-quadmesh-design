include (plugin)

openflipper_plugin(
        DEPS OpenVolumeMesh CoMISo GMM Boost SUITESPARSE
        TYPES TRIANGLEMESH POLYMESH TETRAHEDRALMESH
)