import type { NextPage } from "next";
import React from "react";
import { useRouter } from "next/router";

import * as ch from "@chakra-ui/react";

const Demo: NextPage = () => {
  const [selected, setSelected] = React.useState("none");
  return (
    <>
      <ch.Flex w="100vw" h="full">
        <ch.Flex shadow="md" w="full" h="full" p="10" m="50">
          <ch.VStack w="full">
            <GeometrySelector selected={selected} setSelected={setSelected} />
            <FormNavigation previous="/" next={`/${selected}/`} />
          </ch.VStack>
        </ch.Flex>
      </ch.Flex>
    </>
  );
};

export default Demo;

const GeometrySelector = ({ selected, setSelected }) => {
  const image_links = {
    none: "",
    sphere: "",
    onion:
      "https://user-images.githubusercontent.com/63464503/124515938-880a8f80-ddd8-11eb-9439-409381b5124a.png",
    "sphere-cut":
      "https://user-images.githubusercontent.com/63464503/137699819-3dad4fb6-7e76-4a5c-89b4-86924a62105c.png",
  };

  return (
    <>
      <ch.FormControl>
        <ch.FormLabel>Choose Geometry</ch.FormLabel>
        <ch.Select
          placeholder="Select geometry"
          onChange={(event) => setSelected(event.target.value)}
        >
          <option value="sphere">
            Sphere (image not available) (option 404's)
          </option>
          <option value="onion">Onion</option>
          <option value="sphere-cut">A Cut Sphere (option 404's)</option>
        </ch.Select>

        <ch.FormHelperText>
          In the full version, you programatically define any geometry.
        </ch.FormHelperText>

        <ch.Flex p={10}>
          {
            image_links[selected] ? (
              <ch.Image src={image_links[selected]} alt={selected} />
            ) : (
              <></>
            )
          }
        </ch.Flex>
      </ch.FormControl>
    </>
  );
};

const FormNavigation = ({ previous, next }) => {
  const router = useRouter();
  return (
    <>
      <ch.HStack w="full" align="left">
        <ch.Spacer />
        <ch.Button onClick={() => router.push(previous)}>Previous</ch.Button>
        <ch.Button colorScheme="facebook" onClick={() => router.push(next)}>
          Next
        </ch.Button>
      </ch.HStack>
    </>
  );
};
