
{
  description = "Environment with named dev environments and inherited shell specs.";

  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs/nixos-unstable";
  };

  outputs = { self, nixpkgs }:
  let
    system = "x86_64-linux";
    pkgs = nixpkgs.legacyPackages.${system};
  in
  {
    devShells = {
      ${system} = {
        default = pkgs.mkShell {
          buildInputs = [
            pkgs.zathura
            pkgs.typst
          ];
          shellHook = ''
            # Start zathura
            ${pkgs.zathura}/bin/zathura main.pdf > zathura.log 2>&1 &
            ZATHURA_PID=$!

            ${pkgs.kitty}/bin/kitty ${pkgs.typst}/bin/typst watch main.typ main.pdf 2>&1 &
            TYPST_PID=$!

            # Clean up on shell exit
            trap 'kill "$ZATHURA_PID" "$TYPST_PID"' EXIT

            hx main.typ
          '';
        };
      };
    };
  };
}
