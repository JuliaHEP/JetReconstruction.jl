# A very simple interrogator to quickly get the contents of an EDM4hep ROOT
# file - useful for debugging and developing!
using EDM4hep
using EDM4hep.RootIO
using JetReconstruction

function main()
    input_file = ARGS[1]
    reader = RootIO.Reader(input_file)
    events = RootIO.get(reader, "events")
    evt = events[1]
    # See everything!
    show(IOContext(stdout, :limit => false, :compact => false, :displaysize => (1000, 200)),
         reader)
end

main()
