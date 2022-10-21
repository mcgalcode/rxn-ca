# rxn-ca - A lattice model for simulating solid state reactions

## Debugging

### grpcio

On M1 Macs, importing code related to the reaction-network can cause breakages. The suggested action from the error message, reproduced below, should fix any issues.

```
Failed to import grpc on Apple Silicon. On Apple Silicon machines, try `pip uninstall grpcio; conda install grpcio`. Check out https://docs.ray.io/en/master/ray-overview/installation.html#m1-mac-apple-silicon-support for more details.
```
